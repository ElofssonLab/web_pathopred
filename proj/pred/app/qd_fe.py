#!/usr/bin/env python
# Description: daemon to submit jobs and retrieve results to/from remote
#              servers
# 
import os
import sys
import site

rundir = os.path.dirname(os.path.realpath(__file__))
webserver_root = os.path.realpath("%s/../../../"%(rundir))

activate_env="%s/env/bin/activate_this.py"%(webserver_root)
exec(compile(open(activate_env, "rb").read(), activate_env, 'exec'), dict(__file__=activate_env))
#Add the site-packages of the virtualenv
site.addsitedir("%s/env/lib/python2.7/site-packages/"%(webserver_root))
sys.path.append("%s/env/lib/python2.7/site-packages/"%(webserver_root))
sys.path.append("/usr/local/lib/python2.7/dist-packages")

from libpredweb import myfunc
from libpredweb import dataprocess
from libpredweb import webserver_common as webcom
from libpredweb import qd_fe_common as qdcom

from datetime import datetime
from dateutil import parser as dtparser
from pytz import timezone

import time
import requests
import json
import urllib.request, urllib.parse, urllib.error
import shutil
import hashlib
import subprocess
from suds.client import Client
import numpy

from geoip import geolite2
import pycountry

TZ = webcom.TZ
os.environ['TZ'] = TZ
time.tzset()

vip_user_list = [
        "nanjiang.shu@scilifelab.se"
        ]

# make sure that only one instance of the script is running
# this code is working 
progname = os.path.basename(__file__)
rootname_progname = os.path.splitext(progname)[0]
lockname = os.path.realpath(__file__).replace(" ", "").replace("/", "-")
import fcntl
lock_file = "/tmp/%s.lock"%(lockname)
fp = open(lock_file, 'w')
try:
    fcntl.lockf(fp, fcntl.LOCK_EX | fcntl.LOCK_NB)
except IOError:
    print("Another instance of %s is running"%(progname), file=sys.stderr)
    sys.exit(1)

contact_email = "nanjiang.shu@scilifelab.se"

threshold_logfilesize = 20*1024*1024

usage_short="""
Usage: %s
"""%(sys.argv[0])

usage_ext="""
Description:
    Daemon to submit jobs and retrieve results to/from remote servers
    run periodically
    At the end of each run generate a runlog file with the status of all jobs

OPTIONS:
  -h, --help    Print this help message and exit

Created 2015-03-25, updated 2015-03-25, Nanjiang Shu
"""
usage_exp="""
"""

basedir = os.path.realpath("%s/.."%(rundir)) # path of the application, i.e. pred/
path_static = "%s/static"%(basedir)
path_log = "%s/static/log"%(basedir)
path_stat = "%s/stat"%(path_log)
path_result = "%s/static/result"%(basedir)
path_cache = "%s/static/result/cache"%(basedir)
name_cachedir = 'cache'

# format of the computenodefile is 
# each line is a record and contains two items
# hostname MAX_ALLOWED_PARALLEL_JOBS
computenodefile = "%s/config/computenode.txt"%(basedir)

gen_errfile = "%s/static/log/%s.err"%(basedir, progname)
gen_logfile = "%s/static/log/%s.log"%(basedir, progname)
black_iplist_file = "%s/config/black_iplist.txt"%(basedir)
finished_date_db = "%s/cached_job_finished_date.sqlite3"%(path_log)
vip_email_file = "%s/config/vip_email.txt"%(basedir)

def PrintHelp(fpout=sys.stdout):#{{{
    print(usage_short, file=fpout)
    print(usage_ext, file=fpout)
    print(usage_exp, file=fpout)#}}}

def CreateRunJoblog(path_result, submitjoblogfile, runjoblogfile,#{{{
        finishedjoblogfile, loop):
    myfunc.WriteFile("CreateRunJoblog...\n", gen_logfile, "a", True)
    # Read entries from submitjoblogfile, checking in the result folder and
    # generate two logfiles: 
    #   1. runjoblogfile 
    #   2. finishedjoblogfile
    # when loop == 0, for unfinished jobs, re-generate finished_seqs.txt
    hdl = myfunc.ReadLineByBlock(submitjoblogfile)
    if hdl.failure:
        return 1

    finished_jobid_list = []
    finished_job_dict = {}
    if os.path.exists(finishedjoblogfile):
        finished_job_dict = myfunc.ReadFinishedJobLog(finishedjoblogfile)

    new_finished_list = []  # Finished or Failed
    new_runjob_list = []    # Running
    new_waitjob_list = []    # Queued
    lines = hdl.readlines()
    while lines != None:
        for line in lines:
            strs = line.split("\t")
            if len(strs) < 8:
                continue
            submit_date_str = strs[0]
            jobid = strs[1]
            ip = strs[2]
            numseq_str = strs[3]
            jobname = strs[5]
            email = strs[6].strip()
            method_submission = strs[7]
            start_date_str = ""
            finish_date_str = ""
            rstdir = "%s/%s"%(path_result, jobid)

            numseq = 1
            try:
                numseq = int(numseq_str)
            except:
                pass

            if jobid in finished_job_dict:
                if os.path.exists(rstdir):
                    li = [jobid] + finished_job_dict[jobid]
                    new_finished_list.append(li)
                continue


            status = get_job_status(jobid)

            starttagfile = "%s/%s"%(rstdir, "runjob.start")
            finishtagfile = "%s/%s"%(rstdir, "runjob.finish")
            if os.path.exists(starttagfile):
                start_date_str = myfunc.ReadFile(starttagfile).strip()
            if os.path.exists(finishtagfile):
                finish_date_str = myfunc.ReadFile(finishtagfile).strip()

            li = [jobid, status, jobname, ip, email, numseq_str,
                    method_submission, submit_date_str, start_date_str,
                    finish_date_str]
            if status in ["Finished", "Failed"]:
                new_finished_list.append(li)

            # single-sequence job submitted from the web-page will be
            # submmitted by suq
            UPPER_WAIT_TIME_IN_SEC = 60
            isValidSubmitDate = True
            try:
                submit_date = datetime.datetime.strptime(submit_date_str, "%Y-%m-%d %H:%M:%S")
            except ValueError:
                isValidSubmitDate = False

            if isValidSubmitDate:
                current_time = datetime.datetime.now()
                timeDiff = current_time - submit_date
                queuetime_in_sec = timeDiff.seconds
            else:
                queuetime_in_sec = UPPER_WAIT_TIME_IN_SEC + 1

            if numseq > 1 or method_submission == "wsdl" or queuetime_in_sec > UPPER_WAIT_TIME_IN_SEC:
                if status == "Running":
                    new_runjob_list.append(li)
                elif status == "Wait":
                    new_waitjob_list.append(li)
        lines = hdl.readlines()
    hdl.close()

# re-write logs of finished jobs
    li_str = []
    for li in new_finished_list:
        li_str.append("\t".join(li))
    if len(li_str)>0:
        myfunc.WriteFile("\n".join(li_str)+"\n", finishedjoblogfile, "w", True)
    else:
        myfunc.WriteFile("", finishedjoblogfile, "w", True)
# re-write logs of finished jobs for each IP
    new_finished_dict = {}
    for li in new_finished_list:
        ip = li[3]
        if not ip in new_finished_dict:
            new_finished_dict[ip] = []
        new_finished_dict[ip].append(li)
    for ip in new_finished_dict:
        finished_list_for_this_ip = new_finished_dict[ip]
        divide_finishedjoblogfile = "%s/divided/%s_finished_job.log"%(path_log,
                ip)
        li_str = []
        for li in finished_list_for_this_ip:
            li_str.append("\t".join(li))
        if len(li_str)>0:
            myfunc.WriteFile("\n".join(li_str)+"\n", divide_finishedjoblogfile, "w", True)
        else:
            myfunc.WriteFile("", divide_finishedjoblogfile, "w", True)

# update allfinished jobs
    allfinishedjoblogfile = "%s/all_finished_job.log"%(path_log)
    allfinished_jobid_set = set(myfunc.ReadIDList2(allfinishedjoblogfile, col=0, delim="\t"))
    li_str = []
    for li in new_finished_list:
        jobid = li[0]
        if not jobid in allfinished_jobid_set:
            li_str.append("\t".join(li))
    if len(li_str)>0:
        myfunc.WriteFile("\n".join(li_str)+"\n", allfinishedjoblogfile, "a", True)

# write logs of running and queuing jobs
# the queuing jobs are sorted in descending order by the suq priority
# frist get numseq_this_user for each jobs
# format of numseq_this_user: {'jobid': numseq_this_user}
    numseq_user_dict = GetNumSeqSameUserDict(new_runjob_list + new_waitjob_list)

# now append numseq_this_user and priority score to new_waitjob_list and
# new_runjob_list

    for joblist in [new_waitjob_list, new_runjob_list]:
        for li in joblist:
            jobid = li[0]
            ip = li[3]
            email = li[4].strip()
            rstdir = "%s/%s"%(path_result, jobid)
            outpath_result = "%s/%s"%(rstdir, jobid)

            # if loop == 0 , for new_waitjob_list and new_runjob_list
            # re-generate finished_seqs.txt
            if loop == 0 and os.path.exists(outpath_result):#{{{
                finished_seq_file = "%s/finished_seqs.txt"%(outpath_result)
                finished_idx_file = "%s/finished_seqindex.txt"%(rstdir)
                finished_idx_set = set([])

                finished_seqs_idlist = []
                if os.path.exists(finished_seq_file):
                    finished_seqs_idlist = myfunc.ReadIDList2(finished_seq_file, col=0, delim="\t")
                finished_seqs_idset = set(finished_seqs_idlist)
                finished_info_list = []
                queryfile = "%s/query.fa"%(rstdir)
                (seqidlist, seqannolist, seqlist) = myfunc.ReadFasta(queryfile)
                try:
                    dirlist = os.listdir(outpath_result)
                    for dd in dirlist:
                        if dd.find("seq_") == 0:
                            origIndex_str = dd.split("_")[1]
                            finished_idx_set.add(origIndex_str)

                        if dd.find("seq_") == 0 and dd not in finished_seqs_idset:
                            origIndex = int(dd.split("_")[1])
                            outpath_this_seq = "%s/%s"%(outpath_result, dd)
                            timefile = "%s/time.txt"%(outpath_this_seq)
                            runtime1 = 0.0
                            seq = seqlist[origIndex]
                            description = seqannolist[origIndex]
                            if os.path.exists(timefile):
                                txt = myfunc.ReadFile(timefile).strip()
                                ss2 = txt.split(";")
                                try:
                                    runtime = float(ss2[1])
                                except:
                                    runtime = runtime1
                                    pass
                            else:
                                runtime = runtime1

                            info_finish = [ dd, str(len(seq)), 
                                    None, None,
                                    "newrun", str(runtime), description]
                            finished_info_list.append("\t".join(info_finish))
                except:
                    date_str = time.strftime("%Y-%m-%d %H:%M:%S")
                    myfunc.WriteFile("[Date: %s] Failed to os.listdir(%s)\n"%(date_str, outpath_result), gen_errfile, "a", True)
                    raise
                if len(finished_info_list)>0:
                    myfunc.WriteFile("\n".join(finished_info_list)+"\n", finished_seq_file, "a", True)
                if len(finished_idx_set) > 0:
                    myfunc.WriteFile("\n".join(list(finished_idx_set))+"\n", finished_idx_file, "w", True)
                else:
                    myfunc.WriteFile("", finished_idx_file, "w", True)


            #}}}

            try:
                numseq = int(li[5])
            except:
                numseq = 1
                pass
            try:
                numseq_this_user = numseq_user_dict[jobid]
            except:
                numseq_this_user = numseq
                pass
            # note that the priority is deducted by numseq so that for jobs
            # from the same user, jobs with fewer sequences are placed with
            # higher priority
            priority = myfunc.GetSuqPriority(numseq_this_user) - numseq

            if ip in g_params['blackiplist']:
                priority = priority/1000.0

            if email in vip_user_list:
                numseq_this_user = 1
                priority = 999999999.0
                myfunc.WriteFile("email %s in vip_user_list\n"%(email), gen_logfile, "a", True)

            li.append(numseq_this_user)
            li.append(priority)


    # sort the new_waitjob_list in descending order by priority
    new_waitjob_list = sorted(new_waitjob_list, key=lambda x:x[11], reverse=True)
    new_runjob_list = sorted(new_runjob_list, key=lambda x:x[11], reverse=True)

    # write to runjoblogfile
    li_str = []
    for joblist in [new_waitjob_list, new_runjob_list]:
        for li in joblist:
            li2 = li[:10]+[str(li[10]), str(li[11])]
            li_str.append("\t".join(li2))
#     print "write to", runjoblogfile
#     print "\n".join(li_str)
    if len(li_str)>0:
        myfunc.WriteFile("\n".join(li_str)+"\n", runjoblogfile, "w", True)
    else:
        myfunc.WriteFile("", runjoblogfile, "w", True)

#}}}
def SubmitJob(jobid,cntSubmitJobDict, numseq_this_user):#{{{
# for each job rstdir, keep three log files, 
# 1.seqs finished, finished_seq log keeps all information, finished_index_log
#   can be very compact to speed up reading, e.g.
#   1-5 7-9 etc
# 2.seqs queued remotely , format:
#       index node remote_jobid
# 3. format of the torun_idx_file
#    origIndex

    rmsg = ""
    myfunc.WriteFile("SubmitJob for %s, numseq_this_user=%d\n"%(jobid,
        numseq_this_user), gen_logfile, "a", True)
    rstdir = "%s/%s"%(path_result, jobid)
    outpath_result = "%s/%s"%(rstdir, jobid)
    if not os.path.exists(outpath_result):
        os.mkdir(outpath_result)

    finished_idx_file = "%s/finished_seqindex.txt"%(rstdir)
    failed_idx_file = "%s/failed_seqindex.txt"%(rstdir)
    remotequeue_idx_file = "%s/remotequeue_seqindex.txt"%(rstdir)
    torun_idx_file = "%s/torun_seqindex.txt"%(rstdir) # ordered seq index to run
    cnttry_idx_file = "%s/cntsubmittry_seqindex.txt"%(rstdir)#index file to keep log of tries

    errfile = "%s/%s"%(rstdir, "runjob.err")
    finished_seq_file = "%s/finished_seqs.txt"%(outpath_result)
    tmpdir = "%s/tmpdir"%(rstdir)
    qdinittagfile = "%s/runjob.qdinit"%(rstdir)
    failedtagfile = "%s/%s"%(rstdir, "runjob.failed")
    starttagfile = "%s/%s"%(rstdir, "runjob.start")
    fafile = "%s/query.fa"%(rstdir)
    split_seq_dir = "%s/splitaa"%(tmpdir)
    forceruntagfile = "%s/forcerun"%(rstdir)
    variant_file = "%s/variants.fa"%(rstdir)

    finished_idx_list = []
    failed_idx_list = []    # [origIndex]
    if os.path.exists(finished_idx_file):
        finished_idx_list = list(set(myfunc.ReadIDList(finished_idx_file)))
    if os.path.exists(failed_idx_file):
        failed_idx_list = list(set(myfunc.ReadIDList(failed_idx_file)))

    processed_idx_set = set(finished_idx_list) | set(failed_idx_list)

    jobinfofile = "%s/jobinfo"%(rstdir)
    jobinfo = ""
    if os.path.exists(jobinfofile):
        jobinfo = myfunc.ReadFile(jobinfofile).strip()
    jobinfolist = jobinfo.split("\t")
    email = ""
    if len(jobinfolist) >= 8:
        email = jobinfolist[6]
        method_submission = jobinfolist[7]

    if not os.path.exists(tmpdir):
        os.mkdir(tmpdir)

    # the first time when the this jobid is processed, do the following
    # 1. generate a file with sorted seqindex
    # 2. generate splitted sequence files named by the original seqindex
    if not os.path.exists(qdinittagfile): #initialization#{{{

        init_finished_idx_list = [] # [origIndex]
        if os.path.exists(finished_idx_file):
            init_finished_idx_list = myfunc.ReadIDList(finished_idx_file)

        init_finished_idx_set = set(init_finished_idx_list)
        if g_params['DEBUG']:
            myfunc.WriteFile("DEBUG: init_finished_idx_list = %s"%(str(init_finished_idx_list)), gen_logfile, "a", True)

        # ==== 1.dealing with cached results 
        (seqIDList, seqAnnoList, seqList) = myfunc.ReadFasta(fafile)
        if len(seqIDList) <= 0:
            date_str = time.strftime("%Y-%m-%d %H:%M:%S")
            myfunc.WriteFile(date_str, failedtagfile, "w", True)
            myfunc.WriteFile("Read query seq file failed. Zero sequence read in.\n", errfile, "a", True)
            return 1
        toRunDict = {}
        if os.path.exists(forceruntagfile):
            for i in range(len(seqIDList)):
                toRunDict[i] = [seqList[i], 0, seqAnnoList[i]]
        else:
            for i in range(len(seqIDList)):
                if not str(i) in init_finished_idx_set:
                    toRunDict[i] = [seqList[i], 0, seqAnnoList[i]] #init value for numTM is 0


        sortedlist = sorted(list(toRunDict.items()), key=lambda x:x[1][1], reverse=True)

        # Write splitted fasta file and write a torunlist.txt
        if not os.path.exists(split_seq_dir):
            os.mkdir(split_seq_dir)

        torun_index_str_list = [str(x[0]) for x in sortedlist]
        if len(torun_index_str_list)>0:
            myfunc.WriteFile("\n".join(torun_index_str_list)+"\n", torun_idx_file, "w", True)
        else:
            myfunc.WriteFile("", torun_idx_file, "w", True)

        # write cnttry file for each jobs to run
        cntTryDict = {}
        for idx in torun_index_str_list:
            cntTryDict[int(idx)] = 0
        json.dump(cntTryDict, open(cnttry_idx_file, "w"))

        for item in sortedlist:
            origIndex = item[0]
            seq = item[1][0]
            description = item[1][2]
            seqfile_this_seq = "%s/%s"%(split_seq_dir, "query_%d.fa"%(origIndex))
            seqcontent = ">%s\n%s\n"%(description, seq)
            myfunc.WriteFile(seqcontent, seqfile_this_seq, "w", True)
        # qdinit file is written at the end of initialization, to make sure
        # that initialization is either not started or completed
        date_str = time.strftime("%Y-%m-%d %H:%M:%S")
        myfunc.WriteFile(date_str, qdinittagfile, "w", True)
#}}}


    # 3. try to submit the job 
    if os.path.exists(forceruntagfile):
        isforcerun = "True"
    else:
        isforcerun = "False"
    toRunIndexList = [] # index in str
    processedIndexSet = set([]) #seq index set that are already processed
    submitted_loginfo_list = []
    if os.path.exists(torun_idx_file):
        toRunIndexList = myfunc.ReadIDList(torun_idx_file)
        # unique the list but keep the order
        toRunIndexList = myfunc.uniquelist(toRunIndexList)
    if len(toRunIndexList) > 0:
        iToRun = 0
        numToRun = len(toRunIndexList)
        for node in cntSubmitJobDict:
            if iToRun >= numToRun:
                break
            wsdl_url = "http://%s/pred/api_submitseq/?wsdl"%(node)
            try:
                myclient = Client(wsdl_url, cache=None, timeout=30)
            except:
                date_str = time.strftime("%Y-%m-%d %H:%M:%S")
                myfunc.WriteFile("[Date: %s] Failed to access %s\n"%(date_str,
                    wsdl_url), gen_errfile, "a", True)
                break

            [cnt, maxnum] = cntSubmitJobDict[node]
            cnttry = 0
            while cnt < maxnum and iToRun < numToRun:
                origIndex = int(toRunIndexList[iToRun])
                seqfile_this_seq = "%s/%s"%(split_seq_dir, "query_%d.fa"%(origIndex))
                # ignore already existing query seq, this is an ugly solution,
                # the generation of torunindexlist has a bug
                outpath_this_seq = "%s/%s"%(outpath_result, "seq_%d"%origIndex)
                if os.path.exists(outpath_this_seq):
                    iToRun += 1
                    continue


                if g_params['DEBUG']:
                    myfunc.WriteFile("DEBUG: cnt (%d) < maxnum (%d) "\
                            "and iToRun(%d) < numToRun(%d)"%(cnt, maxnum, iToRun, numToRun), gen_logfile, "a", True)
                fastaseq = ""
                seqid = ""
                seqanno = ""
                seq = ""
                if not os.path.exists(seqfile_this_seq):
                    all_seqfile = "%s/query.fa"%(rstdir)
                    try:
                        (allseqidlist, allannolist, allseqlist) = myfunc.ReadFasta(all_seqfile)
                        seqid = allseqidlist[origIndex]
                        seqanno = allannolist[origIndex]
                        seq = allseqlist[origIndex]
                        fastaseq = ">%s\n%s\n"%(seqanno, seq)
                    except:
                        pass
                else:
                    fastaseq = myfunc.ReadFile(seqfile_this_seq)#seq text in fasta format
                    (seqid, seqanno, seq) = myfunc.ReadSingleFasta(seqfile_this_seq)


                isSubmitSuccess = False
                if len(seq) > 0:
                    query_para = {}
                    query_para['name_software'] = "docker_pathopred"
                    #include variants as a query parameter
                    variant_text = myfunc.ReadFile(variant_file)
                    query_para['variants'] = variant_text
                    #also include the identifier name as a query parameter
                    query_para['identifier_name'] = seqid

                    para_str = json.dumps(query_para, sort_keys=True)
                    jobname = ""
                    if not email in vip_user_list:
                        useemail = ""
                    else:
                        useemail = email
                    try:
                        myfunc.WriteFile("\tSubmitting seq %4d "%(origIndex),
                                gen_logfile, "a", True)
                        rtValue = myclient.service.submitjob_remote(fastaseq, para_str,
                                jobname, useemail, str(numseq_this_user), isforcerun)
                    except:
                        date_str = time.strftime("%Y-%m-%d %H:%M:%S")
                        myfunc.WriteFile("[Date: %s] Failed to run myclient.service.submitjob_remote\n"%(date_str), gen_errfile, "a", True)
                        rtValue = []
                        pass

                    cnttry += 1
                    if len(rtValue) >= 1:

                        strs = rtValue[0]
                        if len(strs) >=5:
                            remote_jobid = strs[0]
                            result_url = strs[1]
                            numseq_str = strs[2]
                            errinfo = strs[3]
                            warninfo = strs[4]
                            if remote_jobid != "None" and remote_jobid != "":
                                isSubmitSuccess = True
                                epochtime = time.time()
                                # 6 fields in the file remotequeue_idx_file
                                txt =  "%d\t%s\t%s\t%s\t%s\t%f"%( origIndex,
                                        node, remote_jobid, seqanno, seq,
                                        epochtime)
                                submitted_loginfo_list.append(txt)
                                cnttry = 0  #reset cnttry to zero
                        else:
                            date_str = time.strftime("%Y-%m-%d %H:%M:%S")
                            myfunc.WriteFile("[Date: %s] bad wsdl return value\n"%(date_str), gen_errfile, "a", True)
                if isSubmitSuccess:
                    cnt += 1
                    myfunc.WriteFile(" succeeded\n", gen_logfile, "a", True)
                else:
                    myfunc.WriteFile(" failed\n", gen_logfile, "a", True)

                if isSubmitSuccess or cnttry >= g_params['MAX_SUBMIT_TRY']:
                    iToRun += 1
                    processedIndexSet.add(str(origIndex))
                    if g_params['DEBUG']:
                        myfunc.WriteFile("DEBUG: jobid %s processedIndexSet.add(str(%d))\n"%(jobid, origIndex), gen_logfile, "a", True)
            # update cntSubmitJobDict for this node
            cntSubmitJobDict[node] = [cnt, maxnum]

    # finally, append submitted_loginfo_list to remotequeue_idx_file 
    if len(submitted_loginfo_list)>0:
        myfunc.WriteFile("\n".join(submitted_loginfo_list)+"\n", remotequeue_idx_file, "a", True)
    # update torun_idx_file
    newToRunIndexList = []
    for idx in toRunIndexList:
        if not idx in processedIndexSet:
            newToRunIndexList.append(idx)
    if g_params['DEBUG']:
        myfunc.WriteFile("DEBUG: jobid %s, newToRunIndexList="%(jobid) + " ".join( newToRunIndexList)+"\n", gen_logfile, "a", True)

    if len(newToRunIndexList)>0:
        myfunc.WriteFile("\n".join(newToRunIndexList)+"\n", torun_idx_file, "w", True)
    else:
        myfunc.WriteFile("", torun_idx_file, "w", True)

    return 0
#}}}
def GetResult(jobid):#{{{
    # retrieving result from the remote server for this job
    myfunc.WriteFile("GetResult for %s.\n" %(jobid), gen_logfile, "a", True)
    rstdir = "%s/%s"%(path_result, jobid)
    outpath_result = "%s/%s"%(rstdir, jobid)
    if not os.path.exists(outpath_result):
        os.mkdir(outpath_result)

    remotequeue_idx_file = "%s/remotequeue_seqindex.txt"%(rstdir)

    torun_idx_file = "%s/torun_seqindex.txt"%(rstdir) # ordered seq index to run
    finished_idx_file = "%s/finished_seqindex.txt"%(rstdir)
    failed_idx_file = "%s/failed_seqindex.txt"%(rstdir)

    starttagfile = "%s/%s"%(rstdir, "runjob.start")
    cnttry_idx_file = "%s/cntsubmittry_seqindex.txt"%(rstdir)#index file to keep log of tries
    tmpdir = "%s/tmpdir"%(rstdir)
    finished_seq_file = "%s/finished_seqs.txt"%(outpath_result)

    if not os.path.exists(tmpdir):
        os.mkdir(tmpdir)

    finished_info_list = [] #[info for finished record]
    finished_idx_list = [] # [origIndex]
    failed_idx_list = []    # [origIndex]
    resubmit_idx_list = []  # [origIndex]
    keep_queueline_list = [] # [line] still in queue

    cntTryDict = {}
    if os.path.exists(cnttry_idx_file):
        with open(cnttry_idx_file, 'r') as fpin:
            cntTryDict = json.load(fpin)

    # in case of missing queries, if remotequeue_idx_file is empty  but the job
    # is still not finished, force re-creating torun_idx_file
    if ((not os.path.exists(remotequeue_idx_file) or
        os.path.getsize(remotequeue_idx_file)<1)):
        idlist1 = []
        idlist2 = []
        if os.path.exists(finished_idx_file):
           idlist1 =  myfunc.ReadIDList(finished_idx_file)
        if os.path.exists(failed_idx_file):
           idlist2 =  myfunc.ReadIDList(failed_idx_file)

        completed_idx_set = set(idlist1 + idlist2)

        jobinfofile = "%s/jobinfo"%(rstdir)
        jobinfo = myfunc.ReadFile(jobinfofile).strip()
        jobinfolist = jobinfo.split("\t")
        if len(jobinfolist) >= 8:
            numseq = int(jobinfolist[3])

        if len(completed_idx_set) < numseq:
            all_idx_list = [str(x) for x in range(numseq)]
            torun_idx_str_list = list(set(all_idx_list)-completed_idx_set)
            for idx in torun_idx_str_list:
                try:
                    cntTryDict[int(idx)] += 1
                except:
                    cntTryDict[int(idx)] = 1
                    pass
            myfunc.WriteFile("\n".join(torun_idx_str_list)+"\n", torun_idx_file, "w", True)

            if g_params['DEBUG']:
                myfunc.WriteFile("recreate torun_idx_file: jobid = %s, numseq=%d, len(completed_idx_set)=%d, len(torun_idx_str_list)=%d\n"%(jobid, numseq, len(completed_idx_set), len(torun_idx_str_list)), gen_logfile, "a", True)
        else:
            myfunc.WriteFile("", torun_idx_file, "w", True)

    text = ""
    if os.path.exists(remotequeue_idx_file):
        text = myfunc.ReadFile(remotequeue_idx_file)
    if text == "":
        return 1
    lines = text.split("\n")

    nodeSet = set([])
    for i in range(len(lines)):
        line = lines[i]
        if not line or line[0] == "#":
            continue
        strs = line.split("\t")
        if len(strs) != 6:
            continue
        node = strs[1]
        nodeSet.add(node)

    myclientDict = {}
    for node in nodeSet:
        wsdl_url = "http://%s/pred/api_submitseq/?wsdl"%(node)
        try:
            myclient = Client(wsdl_url, cache=None, timeout=30)
            myclientDict[node] = myclient
        except:
            date_str = time.strftime("%Y-%m-%d %H:%M:%S")
            myfunc.WriteFile("[Date: %s] Failed to access %s\n"%(date_str, wsdl_url), gen_errfile, "a", True)
            pass


    for i in range(len(lines)):#{{{
        line = lines[i]

        if g_params['DEBUG']:
            myfunc.WriteFile("Process %s\n"%(line), gen_logfile, "a", True)
        if not line or line[0] == "#":
            continue
        strs = line.split("\t")
        if len(strs) != 6:
            continue
        origIndex = int(strs[0])
        node = strs[1]
        remote_jobid = strs[2]
        description = strs[3]
        seq = strs[4]
        submit_time_epoch = float(strs[5])
        subfoldername_this_seq = "seq_%d"%(origIndex)
        outpath_this_seq = "%s/%s"%(outpath_result, "seq_%d"%origIndex)

        try:
            myclient = myclientDict[node]
        except KeyError:
            continue
        try:
            rtValue = myclient.service.checkjob(remote_jobid)
        except:
            date_str = time.strftime("%Y-%m-%d %H:%M:%S")
            myfunc.WriteFile("[Date: %s] Failed to run myclient.service.checkjob(%s)\n"%(date_str, remote_jobid), gen_errfile, "a", True)
            rtValue = []
            pass
        isSuccess = False
        isFinish_remote = False

        if len(rtValue) >= 1:
            ss2 = rtValue[0]
            if len(ss2)>=3:
                status = ss2[0]
                result_url = ss2[1]
                errinfo = ss2[2]

                if errinfo and errinfo.find("does not exist")!=-1:
                    isFinish_remote = True

                if status == "Finished":#{{{
                    isFinish_remote = True
                    outfile_zip = "%s/%s.zip"%(tmpdir, remote_jobid)
                    isRetrieveSuccess = False
                    myfunc.WriteFile("\tFetching result for %s "%(result_url),
                            gen_logfile, "a", True)
                    if myfunc.IsURLExist(result_url,timeout=5):
                        try:
                            urllib.request.urlretrieve (result_url, outfile_zip)
                            isRetrieveSuccess = True
                            myfunc.WriteFile(" succeeded\n", gen_logfile, "a", True)
                        except Exception as e:
                            myfunc.WriteFile(" failed with %s\n"%(str(e)), gen_logfile, "a", True)
                            pass
                    if os.path.exists(outfile_zip) and isRetrieveSuccess:
                        cmd = ["unzip", outfile_zip, "-d", tmpdir]
                        cmdline = " ".join(cmd)
                        try:
                            rmsg = subprocess.check_output(cmd)
                        except subprocess.CalledProcessError as e:
                            date_str = time.strftime("%Y-%m-%d %H:%M:%S")
                            myfunc.WriteFile("[Date: %s] cmdline=%s\nerrmsg=%s\n"%(
                                    date_str, cmdline, str(e)), gen_errfile, "a", True)
                            pass
                        rst_this_seq = "%s/%s/seq_0"%(tmpdir, remote_jobid)

                        if os.path.islink(outpath_this_seq):
                            os.unlink(outpath_this_seq)
                        elif os.path.exists(outpath_this_seq):
                            shutil.rmtree(outpath_this_seq)

                        if os.path.exists(rst_this_seq) and not os.path.exists(outpath_this_seq):
                            cmd = ["mv","-f", rst_this_seq, outpath_this_seq]
                            cmdline = " ".join(cmd)
                            try:
                                rmsg = subprocess.check_output(cmd)
                            except subprocess.CalledProcessError as e:
                                date_str = time.strftime("%Y-%m-%d %H:%M:%S")
                                myfunc.WriteFile( "[Date: %s] cmdline=%s\nerrmsg=%s\n"%(
                                        date_str, cmdline, str(e)), gen_errfile, "a", True)
                                pass

                            checkfile = "%s/output_predictions"%(outpath_this_seq)
                            if os.path.exists(checkfile):
                                isSuccess = True

                            if isSuccess:
                                # create or update the md5 cache
                                md5_key = hashlib.md5(seq).hexdigest()
                                subfoldername = md5_key[:2]
                                md5_subfolder = "%s/%s"%(path_cache, subfoldername)
                                cachedir = "%s/%s/%s"%(path_cache, subfoldername, md5_key)
                                if os.path.exists(cachedir):
                                    try:
                                        shutil.rmtree(cachedir)
                                    except Exception as e:
                                        myfunc.WriteFile("\tFailed to shutil.rmtree(%s) with %s\n"%(cachedir, str(e)), gen_errfile, "a", True)
                                        pass

                                if not os.path.exists(md5_subfolder):
                                    try:
                                        os.makedirs(md5_subfolder)
                                    except:
                                        pass

                                if os.path.exists(md5_subfolder) and not os.path.exists(cachedir):

                                    cmd = ["mv","-f", outpath_this_seq, cachedir]
                                    cmdline = " ".join(cmd)
                                    try:
                                        subprocess.check_output(cmd)
                                        if g_params['DEBUG_CACHE']:
                                            myfunc.WriteFile("\tDEBUG_CACHE: %s\n"%(cmdline), gen_logfile, "a", True)
                                    except CalledProcessError as e:
                                        print(e)
                                        if g_params['DEBUG_CACHE']:
                                            myfunc.WriteFile("\tDEBUG_CACHE: %s\n"%(str(e)), gen_logfile, "a", True)
                                        pass

                                if not os.path.exists(outpath_this_seq) and os.path.exists(cachedir):
                                    rela_path = os.path.relpath(cachedir, outpath_result) #relative path
                                    try:
                                        if g_params['DEBUG_CACHE']:#{{{
                                            myfunc.WriteFile("\tDEBUG_CACHE: chdir(%s)\n"%(outpath_result), 
                                                    gen_logfile, "a", True)
                                            myfunc.WriteFile("\tDEBUG_CACHE: os.symlink(%s, %s)\n"%(rela_path, 
                                                subfoldername_this_seq), gen_logfile, "a", True)#}}}
                                        os.chdir(outpath_result)
                                        os.symlink(rela_path,  subfoldername_this_seq)
                                    except:
                                        if g_params['DEBUG_CACHE']:#{{{
                                            myfunc.WriteFile("\tDEBUG_CACHE: os.symlink(%s, %s) failed\n"%(rela_path, subfoldername_this_seq), gen_errfile, "a", True)#}}}
                                        pass
                                # delete the data on the remote server
                                try:
                                    rtValue2 = myclient.service.deletejob(remote_jobid)
                                    pass
                                except:
                                    date_str = time.strftime("%Y-%m-%d %H:%M:%S")
                                    myfunc.WriteFile( "[Date: %s] Failed to run myclient.service.deletejob(%s)\n"%(date_str, remote_jobid), gen_errfile, "a", True)
                                    rtValue2 = []
                                    pass

                                logmsg = ""
                                if len(rtValue2) >= 1:
                                    ss2 = rtValue2[0]
                                    if len(ss2) >= 2:
                                        status = ss2[0]
                                        errmsg = ss2[1]
                                        if status == "Succeeded":
                                            logmsg = "Successfully deleted data on %s "\
                                                    "for %s"%(node, remote_jobid)
                                        else:
                                            logmsg = "Failed to delete data on %s for "\
                                                    "%s\nError message:\n%s\n"%(node, remote_jobid, errmsg)
                                else:
                                    logmsg = "Failed to call deletejob %s via WSDL on %s\n"%(remote_jobid, node)

                                # delete the zip file
                                os.remove(outfile_zip)
                                shutil.rmtree("%s/%s"%(tmpdir, remote_jobid))


#}}}
                elif status in ["Failed", "None"]:
                    # the job is failed for this sequence, try to re-submit
                    isFinish_remote = True
                    cnttry = 1
                    try:
                        cnttry = cntTryDict[int(origIndex)]
                    except KeyError:
                        cnttry = 1
                        pass
                    if cnttry < g_params['MAX_RESUBMIT']:
                        resubmit_idx_list.append(str(origIndex))
                        cntTryDict[int(origIndex)] = cnttry+1
                    else:
                        failed_idx_list.append(str(origIndex))
                if status != "Wait" and not os.path.exists(starttagfile):
                    date_str = time.strftime("%Y-%m-%d %H:%M:%S")
                    myfunc.WriteFile(date_str, starttagfile, "w", True)

                if g_params['DEBUG_CACHE']:
                    myfunc.WriteFile("\n", gen_logfile, "a", True)

        if isSuccess:#{{{
            time_now = time.time()
            runtime = 5.0
            runtime1 = time_now - submit_time_epoch #in seconds
            timefile = "%s/time.txt"%(outpath_this_seq)
            if os.path.exists(timefile):
                txt = myfunc.ReadFile(timefile).strip()
                ss2 = txt.split(";")
                try:
                    runtime = float(ss2[1])
                except:
                    runtime = runtime1
                    pass
            else:
                runtime = runtime1

            info_finish = [ "seq_%d"%origIndex, str(len(seq)), 
                    'None', 'None',
                    "newrun", str(runtime), 'None']
            finished_info_list.append("\t".join(info_finish))
            finished_idx_list.append(str(origIndex))#}}}

        if not isFinish_remote:
            time_in_remote_queue = time.time() - submit_time_epoch
            # for jobs queued in the remote queue more than one day (but not
            # running) delete it and try to resubmit it. This solved the
            # problem of dead jobs in the remote server due to server
            # rebooting)
            if status != "Running" and time_in_remote_queue > g_params['MAX_TIME_IN_REMOTE_QUEUE']:
                # delete the remote job on the remote server
                try:
                    rtValue2 = myclient.service.deletejob(remote_jobid)
                except Exception as e:
                    date_str = time.strftime("%Y-%m-%d %H:%M:%S")
                    myfunc.WriteFile( "[Date: %s] Failed to run myclient.service.deletejob(%s) on node %s with msg %s\n"%(date_str, remote_jobid, node, str(e)), gen_logfile, "a", True)
                    rtValue2 = []
                    pass
            else:
                keep_queueline_list.append(line)
#}}}
    #Finally, write log files
    finished_idx_list = list(set(finished_idx_list))
    failed_idx_list = list(set(failed_idx_list))
    resubmit_idx_list = list(set(resubmit_idx_list))


    if len(finished_info_list)>0:
        myfunc.WriteFile("\n".join(finished_info_list)+"\n", finished_seq_file, "a", True)
    if len(finished_idx_list)>0:
        myfunc.WriteFile("\n".join(finished_idx_list)+"\n", finished_idx_file, "a", True)
    if len(failed_idx_list)>0:
        myfunc.WriteFile("\n".join(failed_idx_list)+"\n", failed_idx_file, "a", True)
    if len(resubmit_idx_list)>0:
        myfunc.WriteFile("\n".join(resubmit_idx_list)+"\n", torun_idx_file, "a", True)

    if len(keep_queueline_list)>0:
        myfunc.WriteFile("\n".join(keep_queueline_list)+"\n", remotequeue_idx_file, "w", True);
    else:
        myfunc.WriteFile("", remotequeue_idx_file, "w", True);

    with open(cnttry_idx_file, 'w') as fpout:
        json.dump(cntTryDict, fpout)

    return 0
#}}}

def CheckIfJobFinished(jobid, numseq, email):#{{{
    # check if the job is finished and write tagfiles
    myfunc.WriteFile("CheckIfJobFinished for %s.\n" %(jobid), gen_logfile, "a", True)
    rstdir = "%s/%s"%(path_result, jobid)
    tmpdir = "%s/tmpdir"%(rstdir)
    outpath_result = "%s/%s"%(rstdir, jobid)
    errfile = "%s/%s"%(rstdir, "runjob.err")
    logfile = "%s/%s"%(rstdir, "runjob.log")
    finished_idx_file = "%s/finished_seqindex.txt"%(rstdir)
    failed_idx_file = "%s/failed_seqindex.txt"%(rstdir)
    seqfile = "%s/query.fa"%(rstdir)

    base_www_url_file = "%s/static/log/base_www_url.txt"%(basedir)
    base_www_url = ""

    finished_idx_list = []
    failed_idx_list = []
    if os.path.exists(finished_idx_file):
        finished_idx_list = myfunc.ReadIDList(finished_idx_file)
        finished_idx_list = list(set(finished_idx_list))
    if os.path.exists(failed_idx_file):
        failed_idx_list = myfunc.ReadIDList(failed_idx_file)
        failed_idx_list = list(set(failed_idx_list))

    finishtagfile = "%s/%s"%(rstdir, "runjob.finish")
    failedtagfile = "%s/%s"%(rstdir, "runjob.failed")
    starttagfile = "%s/%s"%(rstdir, "runjob.start")

    num_processed = len(finished_idx_list)+len(failed_idx_list)
    finish_status = "" #["success", "failed", "partly_failed"]
    if num_processed >= numseq:# finished
        if len(failed_idx_list) == 0:
            finish_status = "success"
        elif len(failed_idx_list) >= numseq:
            finish_status = "failed"
        else:
            finish_status = "partly_failed"

        if os.path.exists(base_www_url_file):
            base_www_url = myfunc.ReadFile(base_www_url_file).strip()
        if base_www_url == "":
            base_www_url = "http://pathopred.bioinfo.se"

        date_str = time.strftime("%Y-%m-%d %H:%M:%S")
        date_str_epoch = time.time()
        myfunc.WriteFile(date_str, finishtagfile, "w", True)

        # Now write the text output to a single file
        statfile = "%s/%s"%(outpath_result, "stat.txt")
        resultfile_text = "%s/%s"%(outpath_result, "query.result.txt")
        (seqIDList, seqAnnoList, seqList) = myfunc.ReadFasta(seqfile)
        maplist = []
        for i in range(len(seqIDList)):
            maplist.append("%s\t%d\t%s\t%s"%("seq_%d"%i, len(seqList[i]),
                seqAnnoList[i], seqList[i]))
        start_date_str = myfunc.ReadFile(starttagfile).strip()
        start_date_epoch = datetime.datetime.strptime(start_date_str, "%Y-%m-%d %H:%M:%S").strftime('%s')
        all_runtime_in_sec = float(date_str_epoch) - float(start_date_epoch)

        # now making zip instead (for windows users)
        # note that zip rq will zip the real data for symbolic links
        zipfile = "%s.zip"%(jobid)
        zipfile_fullpath = "%s/%s"%(rstdir, zipfile)
        os.chdir(rstdir)
        cmd = ["zip", "-rq", zipfile, jobid]
        try:
            subprocess.check_output(cmd)
        except subprocess.CalledProcessError as e:
            myfunc.WriteFile(str(e)+"\n", errfile, "a", True)
            pass

        if len(failed_idx_list)>0:
            myfunc.WriteFile(date_str, failedtagfile, "w", True)

        if finish_status == "success":
            shutil.rmtree(tmpdir)

        # send the result to email
        if myfunc.IsValidEmailAddress(email):#{{{

            if os.path.exists(errfile):
                err_msg = myfunc.ReadFile(errfile)

            from_email = "info@pathopred.bioinfo.se"
            to_email = email
            subject = "Your result for Pathopred JOBID=%s"%(jobid)
            if finish_status == "success":
                bodytext = """
    Your result is ready at %s/pred/result/%s

    Thanks for using Pathopred

            """%(base_www_url, jobid)
            elif finish_status == "failed":
                bodytext="""
    We are sorry that your job with jobid %s is failed.

    Please contact %s if you have any questions.

    Attached below is the error message:
    %s
                """%(jobid, contact_email, err_msg)
            else:
                bodytext="""
    Your result is ready at %s/pred/result/%s

    We are sorry that Pathopred failed to predict some sequences of your job.

    Please re-submit the queries that have been failed.

    If you have any further questions, please contact %s.

    Attached below is the error message:
    %s
                """%(base_www_url, jobid, contact_email, err_msg)

            myfunc.WriteFile("Sendmail %s -> %s, %s"% (from_email, to_email, subject), logfile, "a", True)
            rtValue = myfunc.Sendmail(from_email, to_email, subject, bodytext)
            if rtValue != 0:
                myfunc.WriteFile("Sendmail to {} failed with status {}".format(to_email,
                    rtValue), errfile, "a", True)

#}}}
#}}}
def DeleteOldResult(path_result, path_log):#{{{
    """
    Delete jobdirs that are finished > MAX_KEEP_DAYS
    """
    finishedjoblogfile = "%s/finished_job.log"%(path_log)
    finished_job_dict = myfunc.ReadFinishedJobLog(finishedjoblogfile)
    for jobid in finished_job_dict:
        li = finished_job_dict[jobid]
        try:
            finish_date_str = li[8]
        except IndexError:
            finish_date_str = ""
            pass
        if finish_date_str != "":
            isValidFinishDate = True
            try:
                finish_date = datetime.datetime.strptime(finish_date_str, "%Y-%m-%d %H:%M:%S")
            except ValueError:
                isValidFinishDate = False

            if isValidFinishDate:
                current_time = datetime.datetime.now()
                timeDiff = current_time - finish_date
                if timeDiff.days > g_params['MAX_KEEP_DAYS']:
                    rstdir = "%s/%s"%(path_result, jobid)
                    date_str = time.strftime("%Y-%m-%d %H:%M:%S")
                    msg = "\tjobid = %s finished %d days ago (>%d days), delete."%(jobid, timeDiff.days, g_params['MAX_KEEP_DAYS'])
                    myfunc.WriteFile("[Date: %s] "%(date_str)+ msg + "\n", gen_logfile, "a", True)
                    shutil.rmtree(rstdir)
#}}}
def RunStatistics(path_result, path_log):#{{{
# 1. calculate average running time, only for those sequences with time.txt
# show also runtime of type and runtime -vs- seqlength
    myfunc.WriteFile("RunStatistics...\n", gen_logfile, "a", True)
    allfinishedjoblogfile = "%s/all_finished_job.log"%(path_log)
    runtimelogfile = "%s/jobruntime.log"%(path_log)
    runtimelogfile_finishedjobid = "%s/jobruntime_finishedjobid.log"%(path_log)
    submitjoblogfile = "%s/submitted_seq.log"%(path_log)
    if not os.path.exists(path_stat):
        os.mkdir(path_stat)

    allfinishedjobidlist = myfunc.ReadIDList2(allfinishedjoblogfile, col=0, delim="\t")
    runtime_finishedjobidlist = myfunc.ReadIDList(runtimelogfile_finishedjobid)
    toana_jobidlist = list(set(allfinishedjobidlist)-set(runtime_finishedjobidlist))

    for jobid in toana_jobidlist:
        runtimeloginfolist = []
        rstdir = "%s/%s"%(path_result, jobid)
        outpath_result = "%s/%s"%(rstdir, jobid)
        finished_seq_file = "%s/finished_seqs.txt"%(outpath_result)
        lines = []
        if os.path.exists(finished_seq_file):
            lines = myfunc.ReadFile(finished_seq_file).split("\n")
        for line in lines:
            strs = line.split("\t")
            if len(strs)>=7:
                str_seqlen = strs[1]
                str_loc_def = strs[2]
                str_loc_def_score = strs[3]
                source = strs[4]
                if source == "newrun":
                    subfolder = strs[0]
                    timefile = "%s/%s/%s"%(outpath_result, subfolder, "time.txt")
                    if os.path.exists(timefile) and os.path.getsize(timefile)>0:
                        txt = myfunc.ReadFile(timefile).strip()
                        try:
                            ss2 = txt.split(";")
                            runtime_str = ss2[1]
                            if len(ss2) >= 3:
                                database_mode = ss2[2]
                            else:
                                 #set default value of database_mode if it is not available in the timefile
                                database_mode = "PRODRES"  
                            runtimeloginfolist.append("\t".join([jobid, subfolder,
                                source, runtime_str, database_mode, str_seqlen,
                                str_loc_def, str_loc_def_score]))
                        except:
                            sys.stderr.write("bad timefile %s\n"%(timefile))

        if len(runtimeloginfolist)>0:
            # items 
            # jobid, seq_no, newrun_or_cached, runtime, mtd_profile, seqlen, numTM, iShasSP
            myfunc.WriteFile("\n".join(runtimeloginfolist)+"\n",runtimelogfile, "a", True)
        myfunc.WriteFile(jobid+"\n", runtimelogfile_finishedjobid, "a", True)

#2. get numseq_in_job vs count_of_jobs, logscale in x-axis
#   get numseq_in_job vs waiting time (time_start - time_submit)
#   get numseq_in_job vs finish time  (time_finish - time_submit)

    allfinished_job_dict = myfunc.ReadFinishedJobLog(allfinishedjoblogfile)
    countjob_country = {} # countjob_country['country'] = [numseq, numjob, ip_set]
    outfile_numseqjob = "%s/numseq_of_job.stat.txt"%(path_stat)
    outfile_numseqjob_web = "%s/numseq_of_job.web.stat.txt"%(path_stat)
    outfile_numseqjob_wsdl = "%s/numseq_of_job.wsdl.stat.txt"%(path_stat)
    countjob_numseq_dict = {} # count the number jobs for each numseq
    countjob_numseq_dict_web = {} # count the number jobs for each numseq submitted via web
    countjob_numseq_dict_wsdl = {} # count the number jobs for each numseq submitted via wsdl

    waittime_numseq_dict = {}
    waittime_numseq_dict_web = {}
    waittime_numseq_dict_wsdl = {}

    finishtime_numseq_dict = {}
    finishtime_numseq_dict_web = {}
    finishtime_numseq_dict_wsdl = {}

    for jobid in allfinished_job_dict:
        li = allfinished_job_dict[jobid]
        numseq = -1
        try:
            numseq = int(li[4])
        except:
            pass
        try:
            method_submission = li[5]
        except:
            method_submission = ""

        ip = ""
        try:
            ip = li[2]
        except:
            pass

        country = "N/A"           # this is slow
        try:
            match = geolite2.lookup(ip)
            country = pycountry.countries.get(alpha_2=match.country).name
        except:
            pass
        if country != "N/A":
            if not country in countjob_country:
                countjob_country[country] = [0,0,set([])] #[numseq, numjob, ip_set] 
            if numseq != -1:
                countjob_country[country][0] += numseq
            countjob_country[country][1] += 1
            countjob_country[country][2].add(ip)


        submit_date_str = li[6]
        start_date_str = li[7]
        finish_date_str = li[8]

        if numseq != -1:
            if not numseq in  countjob_numseq_dict:
                countjob_numseq_dict[numseq] = 0
            countjob_numseq_dict[numseq] += 1
            if method_submission == "web":
                if not numseq in  countjob_numseq_dict_web:
                    countjob_numseq_dict_web[numseq] = 0
                countjob_numseq_dict_web[numseq] += 1
            if method_submission == "wsdl":
                if not numseq in  countjob_numseq_dict_wsdl:
                    countjob_numseq_dict_wsdl[numseq] = 0
                countjob_numseq_dict_wsdl[numseq] += 1

#           # calculate waittime and finishtime
            isValidSubmitDate = True
            isValidStartDate = True
            isValidFinishDate = True
            try:
                submit_date = datetime.datetime.strptime(submit_date_str, "%Y-%m-%d %H:%M:%S")
            except ValueError:
                isValidSubmitDate = False
            try:
                start_date = datetime.datetime.strptime(start_date_str, "%Y-%m-%d %H:%M:%S")
            except ValueError:
                isValidStartDate = False
            try:
                finish_date = datetime.datetime.strptime(finish_date_str, "%Y-%m-%d %H:%M:%S")
            except ValueError:
                isValidFinishDate = False

            if isValidSubmitDate and isValidStartDate:
                waittime_sec = (start_date - submit_date).total_seconds()
                if not numseq in waittime_numseq_dict:
                    waittime_numseq_dict[numseq] = []
                waittime_numseq_dict[numseq].append(waittime_sec)
                if method_submission == "web":
                    if not numseq in waittime_numseq_dict_web:
                        waittime_numseq_dict_web[numseq] = []
                    waittime_numseq_dict_web[numseq].append(waittime_sec)
                if method_submission == "wsdl":
                    if not numseq in waittime_numseq_dict_wsdl:
                        waittime_numseq_dict_wsdl[numseq] = []
                    waittime_numseq_dict_wsdl[numseq].append(waittime_sec)
            if isValidSubmitDate and isValidFinishDate:
                finishtime_sec = (finish_date - submit_date).total_seconds()
                if not numseq in finishtime_numseq_dict:
                    finishtime_numseq_dict[numseq] = []
                finishtime_numseq_dict[numseq].append(finishtime_sec)
                if method_submission == "web":
                    if not numseq in finishtime_numseq_dict_web:
                        finishtime_numseq_dict_web[numseq] = []
                    finishtime_numseq_dict_web[numseq].append(finishtime_sec)
                if method_submission == "wsdl":
                    if not numseq in finishtime_numseq_dict_wsdl:
                        finishtime_numseq_dict_wsdl[numseq] = []
                    finishtime_numseq_dict_wsdl[numseq].append(finishtime_sec)


    # output countjob by country
    outfile_countjob_by_country = "%s/countjob_by_country.txt"%(path_stat)
    # sort by numseq in descending order
    li_countjob = sorted(list(countjob_country.items()), key=lambda x:x[1][0], reverse=True) 
    li_str = []
    li_str.append("#Country\tNumSeq\tNumJob\tNumIP")
    for li in li_countjob:
        li_str.append("%s\t%d\t%d\t%d"%(li[0], li[1][0], li[1][1], len(li[1][2])))
    myfunc.WriteFile("\n".join(li_str)+"\n", outfile_countjob_by_country, "w", True)

    flist = [outfile_numseqjob, outfile_numseqjob_web, outfile_numseqjob_wsdl  ]
    dictlist = [countjob_numseq_dict, countjob_numseq_dict_web, countjob_numseq_dict_wsdl]
    for i in range(len(flist)):
        dt = dictlist[i]
        outfile = flist[i]
        sortedlist = sorted(list(dt.items()), key = lambda x:x[0])
        try:
            fpout = open(outfile,"w")
            for j in range(len(sortedlist)):
                nseq = sortedlist[j][0]
                count = sortedlist[j][1]
                fpout.write("%d\t%d\n"%(nseq,count))
            fpout.close()
            #plot
            cmd = ["%s/app/plot_numseq_of_job.sh"%(basedir), outfile]
            cmdline = " ".join(cmd)
            try:
                rmsg = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as e:
                date_str = time.strftime("%Y-%m-%d %H:%M:%S")
                myfunc.WriteFile("[Date: %s]"%(date_str)+str(e)+"\n", gen_errfile, "a", True)
                myfunc.WriteFile("[Date: %s] cmdline = %s\n"%(date_str,
                    cmdline), gen_errfile, "a", True)
                pass
        except IOError:
            continue
    cmd = ["%s/app/plot_numseq_of_job_mtp.sh"%(basedir), "-web",
            outfile_numseqjob_web, "-wsdl", outfile_numseqjob_wsdl]
    cmdline = " ".join(cmd)
    try:
        rmsg = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        date_str = time.strftime("%Y-%m-%d %H:%M:%S")
        myfunc.WriteFile("[Date: %s]"%(date_str)+str(e)+"\n", gen_errfile, "a", True)
        myfunc.WriteFile("[Date: %s] cmdline = %s\n"%(date_str,
            cmdline), gen_errfile, "a", True)
        pass

# output waittime vs numseq_of_job
# output finishtime vs numseq_of_job
    outfile_waittime_nseq = "%s/waittime_nseq.stat.txt"%(path_stat)
    outfile_waittime_nseq_web = "%s/waittime_nseq_web.stat.txt"%(path_stat)
    outfile_waittime_nseq_wsdl = "%s/waittime_nseq_wsdl.stat.txt"%(path_stat)
    outfile_finishtime_nseq = "%s/finishtime_nseq.stat.txt"%(path_stat)
    outfile_finishtime_nseq_web = "%s/finishtime_nseq_web.stat.txt"%(path_stat)
    outfile_finishtime_nseq_wsdl = "%s/finishtime_nseq_wsdl.stat.txt"%(path_stat)

    outfile_avg_waittime_nseq = "%s/avg_waittime_nseq.stat.txt"%(path_stat)
    outfile_avg_waittime_nseq_web = "%s/avg_waittime_nseq_web.stat.txt"%(path_stat)
    outfile_avg_waittime_nseq_wsdl = "%s/avg_waittime_nseq_wsdl.stat.txt"%(path_stat)
    outfile_avg_finishtime_nseq = "%s/avg_finishtime_nseq.stat.txt"%(path_stat)
    outfile_avg_finishtime_nseq_web = "%s/avg_finishtime_nseq_web.stat.txt"%(path_stat)
    outfile_avg_finishtime_nseq_wsdl = "%s/avg_finishtime_nseq_wsdl.stat.txt"%(path_stat)

    outfile_median_waittime_nseq = "%s/median_waittime_nseq.stat.txt"%(path_stat)
    outfile_median_waittime_nseq_web = "%s/median_waittime_nseq_web.stat.txt"%(path_stat)
    outfile_median_waittime_nseq_wsdl = "%s/median_waittime_nseq_wsdl.stat.txt"%(path_stat)
    outfile_median_finishtime_nseq = "%s/median_finishtime_nseq.stat.txt"%(path_stat)
    outfile_median_finishtime_nseq_web = "%s/median_finishtime_nseq_web.stat.txt"%(path_stat)
    outfile_median_finishtime_nseq_wsdl = "%s/median_finishtime_nseq_wsdl.stat.txt"%(path_stat)

    flist1 = [ outfile_waittime_nseq , outfile_waittime_nseq_web ,
            outfile_waittime_nseq_wsdl , outfile_finishtime_nseq ,
            outfile_finishtime_nseq_web , outfile_finishtime_nseq_wsdl
            ]

    flist2 = [ outfile_avg_waittime_nseq , outfile_avg_waittime_nseq_web ,
            outfile_avg_waittime_nseq_wsdl , outfile_avg_finishtime_nseq ,
            outfile_avg_finishtime_nseq_web , outfile_avg_finishtime_nseq_wsdl
            ]
    flist3 = [ outfile_median_waittime_nseq , outfile_median_waittime_nseq_web ,
            outfile_median_waittime_nseq_wsdl , outfile_median_finishtime_nseq ,
            outfile_median_finishtime_nseq_web , outfile_median_finishtime_nseq_wsdl
            ]

    dict_list = [
            waittime_numseq_dict , waittime_numseq_dict_web , waittime_numseq_dict_wsdl , finishtime_numseq_dict , finishtime_numseq_dict_web , finishtime_numseq_dict_wsdl
            ]

    for i in range(len(flist1)):
        dt = dict_list[i]
        outfile1 = flist1[i]
        outfile2 = flist2[i]
        outfile3 = flist3[i]
        sortedlist = sorted(list(dt.items()), key = lambda x:x[0])
        try:
            fpout = open(outfile1,"w")
            for j in range(len(sortedlist)):
                nseq = sortedlist[j][0]
                li_time = sortedlist[j][1]
                for k in range(len(li_time)):
                    fpout.write("%d\t%f\n"%(nseq,li_time[k]))
            fpout.close()
        except IOError:
            pass
        try:
            fpout = open(outfile2,"w")
            for j in range(len(sortedlist)):
                nseq = sortedlist[j][0]
                li_time = sortedlist[j][1]
                avg_time = myfunc.FloatDivision(sum(li_time), len(li_time))
                fpout.write("%d\t%f\n"%(nseq,avg_time))
            fpout.close()
        except IOError:
            pass
        try:
            fpout = open(outfile3,"w")
            for j in range(len(sortedlist)):
                nseq = sortedlist[j][0]
                li_time = sortedlist[j][1]
                median_time = numpy.median(li_time)
                fpout.write("%d\t%f\n"%(nseq,median_time))
            fpout.close()
        except IOError:
            pass

    # plotting 
    flist = flist1
    for i in range(len(flist)):
        outfile = flist[i]
        if os.path.exists(outfile):
            cmd = ["%s/app/plot_nseq_waitfinishtime.sh"%(basedir), outfile]
            cmdline = " ".join(cmd)
            try:
                rmsg = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as e:
                date_str = time.strftime("%Y-%m-%d %H:%M:%S")
                myfunc.WriteFile("[Date: %s]"%(date_str)+str(e)+"\n", gen_errfile, "a", True)
                myfunc.WriteFile("[Date: %s] cmdline = %s\n"%(date_str,
                    cmdline), gen_errfile, "a", True)
                pass
    flist = flist2+flist3
    for i in range(len(flist)):
        outfile = flist[i]
        if os.path.exists(outfile):
            cmd = ["%s/app/plot_avg_waitfinishtime.sh"%(basedir), outfile]
            cmdline = " ".join(cmd)
            try:
                rmsg = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as e:
                date_str = time.strftime("%Y-%m-%d %H:%M:%S")
                myfunc.WriteFile("[Date: %s]"%(date_str)+str(e)+"\n", gen_errfile, "a", True)
                myfunc.WriteFile("[Date: %s] cmdline = %s\n"%(date_str,
                    cmdline), gen_errfile, "a", True)
                pass

# get longest predicted seq
# get query with most TM helics
# get query takes the longest time
    extreme_runtimelogfile = "%s/stat/extreme_jobruntime.log"%(path_log)

    longestlength = -1
    mostTM = -1
    longestruntime = -1.0
    line_longestlength = ""
    line_mostTM = ""
    line_longestruntime = ""

#3. get running time vs sequence length
    cntseq = 0
    cnt_hasSP = 0
    outfile_runtime = "%s/length_runtime.stat.txt"%(path_stat)
    outfile_runtime_pfam = "%s/length_runtime.pfam.stat.txt"%(path_stat)
    outfile_runtime_cdd = "%s/length_runtime.cdd.stat.txt"%(path_stat)
    outfile_runtime_uniref = "%s/length_runtime.uniref.stat.txt"%(path_stat)
    outfile_runtime_avg = "%s/length_runtime.stat.avg.txt"%(path_stat)
    outfile_runtime_pfam_avg = "%s/length_runtime.pfam.stat.avg.txt"%(path_stat)
    outfile_runtime_cdd_avg = "%s/length_runtime.cdd.stat.avg.txt"%(path_stat)
    outfile_runtime_uniref_avg = "%s/length_runtime.uniref.stat.avg.txt"%(path_stat)
    li_length_runtime = []
    li_length_runtime_pfam = []
    li_length_runtime_cdd = []
    li_length_runtime_uniref = []
    dict_length_runtime = {}
    dict_length_runtime_pfam = {}
    dict_length_runtime_cdd = {}
    dict_length_runtime_uniref = {}
    li_length_runtime_avg = []
    li_length_runtime_pfam_avg = []
    li_length_runtime_cdd_avg = []
    li_length_runtime_uniref_avg = []

    if os.path.exists(runtimelogfile):
        hdl = myfunc.ReadLineByBlock(runtimelogfile)
        if not hdl.failure:
            lines = hdl.readlines()
            while lines != None:
                for line in lines:
                    strs = line.split("\t")
                    if len(strs) < 8:
                        continue
                    jobid = strs[0]
                    seqidx = strs[1]
                    runtime = -1.0
                    try:
                        runtime = float(strs[3])
                    except:
                        pass
                    mtd_profile = strs[4]
                    lengthseq = -1
                    try:
                        lengthseq = int(strs[5])
                    except:
                        pass

                    numTM = -1
                    try:
                        numTM = int(strs[6])
                    except:
                        pass
                    isHasSP = strs[7]

                    cntseq += 1
                    if isHasSP == "True":
                        cnt_hasSP += 1

                    if runtime > longestruntime:
                        line_longestruntime = line
                        longestruntime = runtime
                    if lengthseq > longestlength:
                        line_longestseq = line
                        longestlength = lengthseq
                    if numTM > mostTM:
                        mostTM = numTM
                        line_mostTM = line

                    if lengthseq != -1:
                        li_length_runtime.append([lengthseq, runtime])
                        if lengthseq not in dict_length_runtime:
                            dict_length_runtime[lengthseq] = []
                        dict_length_runtime[lengthseq].append(runtime)
                        if mtd_profile == "pfam":
                            li_length_runtime_pfam.append([lengthseq, runtime])
                            if lengthseq not in dict_length_runtime_pfam:
                                dict_length_runtime_pfam[lengthseq] = []
                            dict_length_runtime_pfam[lengthseq].append(runtime)
                        elif mtd_profile == "cdd":
                            li_length_runtime_cdd.append([lengthseq, runtime])
                            if lengthseq not in dict_length_runtime_cdd:
                                dict_length_runtime_cdd[lengthseq] = []
                            dict_length_runtime_cdd[lengthseq].append(runtime)
                        elif mtd_profile == "uniref":
                            li_length_runtime_uniref.append([lengthseq, runtime])
                            if lengthseq not in dict_length_runtime_uniref:
                                dict_length_runtime_uniref[lengthseq] = []
                            dict_length_runtime_uniref[lengthseq].append(runtime)
                lines = hdl.readlines()
            hdl.close()

        li_content = []
        try:
            for line in [line_mostTM, line_longestseq, line_longestruntime]:
                li_content.append(line)
            myfunc.WriteFile("\n".join(li_content)+"\n", extreme_runtimelogfile, "w", True)
        except:
            pass

    # get lengthseq -vs- average_runtime
    dict_list = [dict_length_runtime, dict_length_runtime_pfam, dict_length_runtime_cdd, dict_length_runtime_uniref]
    li_list = [li_length_runtime_avg, li_length_runtime_pfam_avg, li_length_runtime_cdd_avg, li_length_runtime_uniref_avg]
    li_sum_runtime = [0.0]*len(dict_list)
    for i in range(len(dict_list)):
        dt = dict_list[i]
        li = li_list[i]
        for lengthseq in dt:
            avg_runtime = sum(dt[lengthseq])/float(len(dt[lengthseq]))
            li.append([lengthseq, avg_runtime])
            li_sum_runtime[i] += sum(dt[lengthseq])

    avg_runtime = myfunc.FloatDivision(li_sum_runtime[0], len(li_length_runtime))
    avg_runtime_pfam = myfunc.FloatDivision(li_sum_runtime[1], len(li_length_runtime_pfam))
    avg_runtime_cdd = myfunc.FloatDivision(li_sum_runtime[2], len(li_length_runtime_cdd))
    avg_runtime_uniref = myfunc.FloatDivision(li_sum_runtime[3], len(li_length_runtime_uniref))

    li_list = [li_length_runtime, li_length_runtime_pfam,
            li_length_runtime_cdd, li_length_runtime_uniref,
            li_length_runtime_avg, li_length_runtime_pfam_avg,
            li_length_runtime_cdd_avg, li_length_runtime_uniref_avg]
    flist = [outfile_runtime, outfile_runtime_pfam, outfile_runtime_cdd,
            outfile_runtime_uniref, outfile_runtime_avg,
            outfile_runtime_pfam_avg, outfile_runtime_cdd_avg,
            outfile_runtime_uniref_avg]
    for i in range(len(flist)):
        outfile = flist[i]
        li = li_list[i]
        sortedlist = sorted(li, key=lambda x:x[0])
        try:
            fpout = open(outfile,"w")
            for j in range(len(sortedlist)):
                lengthseq = sortedlist[j][0]
                runtime = sortedlist[j][1]
                fpout.write("%d\t%f\n"%(lengthseq,runtime))
            fpout.close()
        except IOError:
            continue

    outfile_avg_runtime = "%s/avg_runtime.stat.txt"%(path_stat)
    try:
        fpout = open(outfile_avg_runtime,"w")
        fpout.write("%s\t%f\n"%("All",avg_runtime))
        fpout.write("%s\t%f\n"%("Pfam",avg_runtime_pfam))
        fpout.write("%s\t%f\n"%("CDD",avg_runtime_cdd))
        fpout.write("%s\t%f\n"%("Uniref",avg_runtime_uniref))
        fpout.close()
    except IOError:
        pass
    if os.path.exists(outfile_avg_runtime):
        cmd = ["%s/app/plot_avg_runtime.sh"%(basedir), outfile_avg_runtime]
        cmdline = " ".join(cmd)
        try:
            rmsg = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            date_str = time.strftime("%Y-%m-%d %H:%M:%S")
            myfunc.WriteFile("[Date: %s]"%(date_str)+str(e)+"\n", gen_errfile, "a", True)
            myfunc.WriteFile("[Date: %s] cmdline = %s\n"%(date_str,
                cmdline), gen_errfile, "a", True)
            pass

    flist = [outfile_runtime, outfile_runtime_pfam, outfile_runtime_cdd,
            outfile_runtime_uniref]
    for outfile in flist:
        if os.path.exists(outfile):
            cmd = ["%s/app/plot_length_runtime.sh"%(basedir), outfile]
            cmdline = " ".join(cmd)
            try:
                rmsg = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as e:
                date_str = time.strftime("%Y-%m-%d %H:%M:%S")
                myfunc.WriteFile("[Date: %s]"%(date_str)+str(e)+"\n", gen_errfile, "a", True)
                myfunc.WriteFile("[Date: %s] cmdline = %s\n"%(date_str,
                    cmdline), gen_errfile, "a", True)
                pass

    cmd = ["%s/app/plot_length_runtime_mtp.sh"%(basedir), "-pfam",
            outfile_runtime_pfam, "-cdd", outfile_runtime_cdd, "-uniref",
            outfile_runtime_uniref, "-sep-avg"]
    cmdline = " ".join(cmd)
    try:
        rmsg = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        date_str = time.strftime("%Y-%m-%d %H:%M:%S")
        myfunc.WriteFile("[Date: %s]"%(date_str)+str(e)+"\n", gen_errfile, "a", True)
        myfunc.WriteFile("[Date: %s] cmdline = %s\n"%(date_str,
            cmdline), gen_errfile, "a", True)
        pass


#5. output num-submission time series with different bins (day, week, month, year)
    hdl = myfunc.ReadLineByBlock(submitjoblogfile)
    dict_submit_day = {}  #["name" numjob, numseq, numjob_web, numseq_web,numjob_wsdl, numseq_wsdl]
    dict_submit_week = {}
    dict_submit_month = {}
    dict_submit_year = {}
    if not hdl.failure:
        lines = hdl.readlines()
        while lines != None:
            for line in lines:
                strs = line.split("\t")
                if len(strs) < 8:
                    continue
                submit_date_str = strs[0]
                numseq = 0
                try:
                    numseq = int(strs[3])
                except:
                    pass
                method_submission = strs[7]
                isValidSubmitDate = True
                try:
                    submit_date = datetime.datetime.strptime(submit_date_str, "%Y-%m-%d %H:%M:%S")
                except ValueError:
                    isValidSubmitDate = False
                if isValidSubmitDate:#{{{
                    day_str = submit_date_str.split()[0]
                    (beginning_of_week, end_of_week) = myfunc.week_beg_end(submit_date)
                    week_str = beginning_of_week.strftime("%Y-%m-%d")
                    month_str = submit_date.strftime("%Y-%b")
                    year_str = submit_date.year
                    day = int(day_str.replace("-", ""))
                    week = int(submit_date.strftime("%Y%V"))
                    month = int(submit_date.strftime("%Y%m"))
                    year = int(year_str)
                    if not day in dict_submit_day:
                                                #all   web  wsdl
                        dict_submit_day[day] = [day_str, 0,0,0,0,0,0]
                    if not week in dict_submit_week:
                        dict_submit_week[week] = [week_str, 0,0,0,0,0,0]
                    if not month in dict_submit_month:
                        dict_submit_month[month] = [month_str, 0,0,0,0,0,0]
                    if not year in dict_submit_year:
                        dict_submit_year[year] = [year_str, 0,0,0,0,0,0]
                    dict_submit_day[day][1] += 1
                    dict_submit_day[day][2] += numseq
                    dict_submit_week[week][1] += 1
                    dict_submit_week[week][2] += numseq
                    dict_submit_month[month][1] += 1
                    dict_submit_month[month][2] += numseq
                    dict_submit_year[year][1] += 1
                    dict_submit_year[year][2] += numseq
                    if method_submission == "web":
                        dict_submit_day[day][3] += 1
                        dict_submit_day[day][4] += numseq
                        dict_submit_week[week][3] += 1
                        dict_submit_week[week][4] += numseq
                        dict_submit_month[month][3] += 1
                        dict_submit_month[month][4] += numseq
                        dict_submit_year[year][3] += 1
                        dict_submit_year[year][4] += numseq
                    if method_submission == "wsdl":
                        dict_submit_day[day][5] += 1
                        dict_submit_day[day][6] += numseq
                        dict_submit_week[week][5] += 1
                        dict_submit_week[week][6] += numseq
                        dict_submit_month[month][5] += 1
                        dict_submit_month[month][6] += numseq
                        dict_submit_year[year][5] += 1
                        dict_submit_year[year][6] += numseq
#}}}
            lines = hdl.readlines()
        hdl.close()

    li_submit_day = []
    li_submit_week = []
    li_submit_month = []
    li_submit_year = []
    li_submit_day_web = []
    li_submit_week_web = []
    li_submit_month_web = []
    li_submit_year_web = []
    li_submit_day_wsdl = []
    li_submit_week_wsdl = []
    li_submit_month_wsdl = []
    li_submit_year_wsdl = []
    dict_list = [dict_submit_day, dict_submit_week, dict_submit_month, dict_submit_year]
    li_list = [ li_submit_day, li_submit_week, li_submit_month, li_submit_year,
            li_submit_day_web, li_submit_week_web, li_submit_month_web, li_submit_year_web,
            li_submit_day_wsdl, li_submit_week_wsdl, li_submit_month_wsdl, li_submit_year_wsdl
            ]

    for i in range(len(dict_list)):
        dt = dict_list[i]
        sortedlist = sorted(list(dt.items()), key = lambda x:x[0])
        for j in range(3):
            li = li_list[j*4+i]
            k1 = j*2 +1
            k2 = j*2 +2
            for kk in range(len(sortedlist)):
                items = sortedlist[kk]
                if items[1][k1] > 0 or items[1][k2] > 0:
                    li.append([items[1][0], items[1][k1], items[1][k2]])

    outfile_submit_day = "%s/submit_day.stat.txt"%(path_stat)
    outfile_submit_week = "%s/submit_week.stat.txt"%(path_stat)
    outfile_submit_month = "%s/submit_month.stat.txt"%(path_stat)
    outfile_submit_year = "%s/submit_year.stat.txt"%(path_stat)
    outfile_submit_day_web = "%s/submit_day_web.stat.txt"%(path_stat)
    outfile_submit_week_web = "%s/submit_week_web.stat.txt"%(path_stat)
    outfile_submit_month_web = "%s/submit_month_web.stat.txt"%(path_stat)
    outfile_submit_year_web = "%s/submit_year_web.stat.txt"%(path_stat)
    outfile_submit_day_wsdl = "%s/submit_day_wsdl.stat.txt"%(path_stat)
    outfile_submit_week_wsdl = "%s/submit_week_wsdl.stat.txt"%(path_stat)
    outfile_submit_month_wsdl = "%s/submit_month_wsdl.stat.txt"%(path_stat)
    outfile_submit_year_wsdl = "%s/submit_year_wsdl.stat.txt"%(path_stat)
    flist = [ 
            outfile_submit_day , outfile_submit_week , outfile_submit_month , outfile_submit_year ,
            outfile_submit_day_web , outfile_submit_week_web , outfile_submit_month_web , outfile_submit_year_web ,
            outfile_submit_day_wsdl , outfile_submit_week_wsdl , outfile_submit_month_wsdl , outfile_submit_year_wsdl 
            ]
    for i in range(len(flist)):
        outfile = flist[i]
        li = li_list[i]
        try:
            fpout = open(outfile,"w")
            for j in range(len(li)):     # name    njob   nseq
                fpout.write("%s\t%d\t%d\n"%(li[j][0], li[j][1], li[j][2]))
            fpout.close()
        except IOError:
            pass
        #plot
        if os.path.exists(outfile):
            cmd = ["%s/app/plot_numsubmit.sh"%(basedir), outfile]
            cmdline = " ".join(cmd)
            try:
                rmsg = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as e:
                date_str = time.strftime("%Y-%m-%d %H:%M:%S")
                myfunc.WriteFile("[Date: %s]"%(date_str)+str(e)+"\n", gen_errfile, "a", True)
                myfunc.WriteFile("[Date: %s] cmdline = %s\n"%(date_str,
                    cmdline), gen_errfile, "a", True)
                pass

#}}}

def main(g_params):#{{{
    submitjoblogfile = "%s/submitted_seq.log"%(path_log)
    runjoblogfile = "%s/runjob_log.log"%(path_log)
    finishedjoblogfile = "%s/finished_job.log"%(path_log)

    if not os.path.exists(path_cache):
        os.mkdir(path_cache)

    loop = 0
    while 1:

        if os.path.exists("%s/CACHE_CLEANING_IN_PROGRESS"%(path_result)):#pause when cache cleaning is in progress
            continue

        # load the config file if exists
        configfile = "%s/config/config.json"%(basedir)
        config = {}
        if os.path.exists(configfile):
            text = myfunc.ReadFile(configfile)
            config = json.loads(text)

        if rootname_progname in config:
            g_params.update(config[rootname_progname])

        if os.path.exists(black_iplist_file):
            g_params['blackiplist'] = myfunc.ReadIDList(black_iplist_file)

        avail_computenode = webcom.ReadComputeNode(computenodefile) # return value is a dict
        g_params['vip_user_list'] = myfunc.ReadIDList2(vip_email_file,  col=0)
        num_avail_node = len(avail_computenode)

        webcom.loginfo("loop %d"%(loop), gen_logfile)

        isOldRstdirDeleted = False
        if loop % g_params['STATUS_UPDATE_FREQUENCY'][0] == g_params['STATUS_UPDATE_FREQUENCY'][1]:
            qdcom.RunStatistics_basic(webserver_root, gen_logfile, gen_errfile)
            isOldRstdirDeleted = webcom.DeleteOldResult(path_result, path_log,
                    gen_logfile, MAX_KEEP_DAYS=g_params['MAX_KEEP_DAYS'])
            webcom.CleanServerFile(path_static, gen_logfile, gen_errfile)
            webcom.CleanCachedResult(path_static, name_cachedir,  gen_logfile, gen_errfile)

        webcom.ArchiveLogFile(path_log, threshold_logfilesize=threshold_logfilesize) 

        qdcom.CreateRunJoblog(loop, isOldRstdirDeleted, g_params)

        # Get number of jobs submitted to the remote server based on the
        # runjoblogfile
        runjobidlist = myfunc.ReadIDList2(runjoblogfile,0)
        remotequeueDict = {}
        for node in avail_computenode:
            remotequeueDict[node] = []
        for jobid in runjobidlist:
            rstdir = "%s/%s"%(path_result, jobid)
            remotequeue_idx_file = "%s/remotequeue_seqindex.txt"%(rstdir)
            if os.path.exists(remotequeue_idx_file):
                content = myfunc.ReadFile(remotequeue_idx_file)
                lines = content.split('\n')
                for line in lines:
                    strs = line.split('\t')
                    if len(strs)>=5:
                        node = strs[1]
                        remotejobid = strs[2]
                        if node in remotequeueDict:
                            remotequeueDict[node].append(remotejobid)

        cntSubmitJobDict = webcom.InitCounterSubmitJobDict(avail_computenode, remotequeueDict, g_params['MAX_SUBMIT_JOB_PER_NODE'])

# entries in runjoblogfile includes jobs in queue or running
        hdl = myfunc.ReadLineByBlock(runjoblogfile)
        if not hdl.failure:
            lines = hdl.readlines()
            while lines != None:
                for line in lines:
                    strs = line.split("\t")
                    if len(strs) >= 11:
                        jobid = strs[0]
                        email = strs[4]
                        try:
                            numseq = int(strs[5])
                        except:
                            numseq = 1
                            pass
                        try:
                            numseq_this_user = int(strs[10])
                        except:
                            numseq_this_user = 1
                            pass
                        rstdir = "%s/%s"%(path_result, jobid)
                        finishtagfile = "%s/%s"%(rstdir, "runjob.finish")
                        status = strs[1]


                        webcom.loginfo("CompNodeStatus: %s"%(str(cntSubmitJobDict)), gen_logfile)

                        runjob_lockfile = "%s/%s/%s.lock"%(path_result, jobid, "runjob.lock")
                        if os.path.exists(runjob_lockfile):
                            msg = "runjob_lockfile %s exists, ignore the job %s" %(runjob_lockfile, jobid)
                            webcom.loginfo(msg, gen_logfile)
                            continue

                        if webcom.IsHaveAvailNode(cntSubmitJobDict):
                            if not g_params['DEBUG_NO_SUBMIT']:
                                qdcom.SubmitJob(jobid, cntSubmitJobDict, numseq_this_user, g_params)
                        qdcom.GetResult(jobid, g_params) # the start tagfile is written when got the first result
                        qdcom.CheckIfJobFinished(jobid, numseq, email, g_params)

                lines = hdl.readlines()
            hdl.close()

        webcom.loginfo("sleep for %d seconds"%(g_params['SLEEP_INTERVAL']), gen_logfile)
        time.sleep(g_params['SLEEP_INTERVAL'])
        loop += 1

    return 0
#}}}


def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    g_params['blackiplist'] = []
    g_params['DEBUG'] = False
    g_params['DEBUG_NO_SUBMIT'] = False
    g_params['DEBUG_CACHE'] = False
    g_params['SLEEP_INTERVAL'] = 5    # sleep interval in seconds
    g_params['MAX_SUBMIT_JOB_PER_NODE'] = 200
    g_params['MAX_KEEP_DAYS'] = 60
    g_params['MAX_RESUBMIT'] = 2
    g_params['MAX_SUBMIT_TRY'] = 3
    g_params['MAX_TIME_IN_REMOTE_QUEUE'] = 3600*24 # one day in seconds
    g_params['STATUS_UPDATE_FREQUENCY'] = [500, 50]  # updated by if loop%$1 == $2
    g_params['FORMAT_DATETIME'] = webcom.FORMAT_DATETIME
    g_params['UPPER_WAIT_TIME_IN_SEC'] = 60
    g_params['name_server'] = "PathoPred"
    g_params['path_static'] = path_static
    g_params['path_result'] = path_result
    g_params['path_log'] = path_log
    g_params['path_cache'] = path_cache
    g_params['vip_email_file'] = vip_email_file
    g_params['gen_logfile'] = gen_logfile
    g_params['finished_date_db'] = finished_date_db
    g_params['gen_errfile'] = gen_errfile
    g_params['contact_email'] = contact_email
    g_params['name_cachedir'] = name_cachedir
    g_params['TZ'] = TZ
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    date_str = time.strftime(g_params['FORMAT_DATETIME'])
    print("\n#%s#\n[Date: %s] qd_fe.py restarted"%('='*80,date_str))
    sys.stdout.flush()
    sys.exit(main(g_params))
