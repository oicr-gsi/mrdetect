import os
import sys
import numpy
import gzip

BLACKLIST_PATH="./blacklist.txt.gz"

def build_R1R2_blacklist(SVM_ROC_lines):
    R1R2_dict={}
    for lines in SVM_ROC_lines:
        items=lines.split()
        chrom=items[0]
        loc=items[1]
        read_id=items[4]
        R1R2_key=chrom+"-"+loc+"-"+read_id
        if R1R2_key not in R1R2_dict:
            R1R2_dict[R1R2_key]=[]
        R1R2_dict[R1R2_key].append(items[3])

    R1R2_fail=[]
    for each_R1R2 in R1R2_dict:
        if len(R1R2_dict[each_R1R2])==2:
            if R1R2_dict[each_R1R2][0]!=R1R2_dict[each_R1R2][1]:
                R1R2_fail.append(each_R1R2)

    return R1R2_fail

def run_R1R2(prefilter_reads,R1R2_blacklist):
    
    pass_filter=[]
    for lines in prefilter_reads:
        items=lines.split()
        chrom=items[0]
        loc=items[1]
        read_id=items[4]
        R1R2_key=chrom+"-"+loc+"-"+read_id
        if R1R2_key not in R1R2_blacklist:
            pass_filter.append(lines)

    return pass_filter



def build_blacklist(relevant_sites):

    relevant_list=[]
    for lines in relevant_sites:
        items=lines.split()
        chrom=items[0]
        loc=items[1]
        site_key=chrom+"-"+loc
        relevant_list.append(site_key)
    
    relevant_list=set(relevant_list)
    blacklist_sites=[]
    for lines in gzip.open(BLACKLIST_PATH):
        if lines[0]=="#":
            continue
        items=lines.split()
        chrom=items[0]
        loc=str(int(items[1])-1)
        site_key=chrom+"-"+loc
        if site_key in relevant_list:
            blacklist_sites.append(site_key)

    return blacklist_sites


def run_blacklist(blacklist,relevant_sites):
    
    passed_filter=[]
    for lines in relevant_sites:
        items=lines.split()
        chrom=items[0]
        loc=items[1]
        site_key=chrom+"-"+loc
        if site_key not in blacklist:
            passed_filter.append(lines)
    return passed_filter



def generate_tumor_dict(tumor_VCF):
    tumor_calls=open(tumor_VCF).readlines()

    tumor_calls_dict={}
    for each_site in tumor_calls:
        if each_site[0]=="#":
            continue
        items=each_site.split()
        chrom=items[0]
        loc=items[1]
        ref=items[3]
        alt=items[4]
        location_key=chrom+"-"+loc
        tumor_calls_dict[location_key]=[[ref,alt]]
    
    return tumor_calls_dict


def make_detections(plasma_reads,tumor_sites_dict):
    detected_reads=[]
    for lines in plasma_reads:
        items=lines.split()
        chrom=items[0]
        location=str(int(items[1])+1)
        site_key=chrom+"-"+location
        ref=items[2]
        alt=items[3]
        tumor_ref=tumor_sites_dict[site_key][0][0]
        tumor_alt=tumor_sites_dict[site_key][0][1]
        if ref==tumor_ref and alt==tumor_alt:
            fragment_info=",".join(items)
            detected_reads.append(lines)
            print(lines)
    return detected_reads


def count_sites(reads_list):
    locations=[]
    for lines in reads_list:
        items=lines.split()
        chrom=items[0]
        location=items[1]
        site_key=chrom+"-"+location
        locations.append(site_key)
    return len(numpy.unique(locations))

def count_fragments(reads_list):
    read_ids=[]
    for lines in reads_list:
        items=lines.split()
        chrom=items[0]
        loc=items[1]
        read_id=items[4]
        read_ids.append(chrom+"-"+loc+"-"+read_id)
    return len(numpy.unique(read_ids))

def SVM_filter(unfiltered_reads,threshold):
    passed_filter=[]
    for lines in unfiltered_reads: 
        items=lines.split()
        SVM_score=float(items[12])
        if SVM_score>threshold:
            passed_filter.append(lines)
    return passed_filter

def SVMroc_filter(unfiltered_reads):
    passed_filter=[]
    for lines in unfiltered_reads:
        items=lines.split()
        SVM_pass=items[13]
        if SVM_pass=="1":
            passed_filter.append(lines)
    return passed_filter    

def main():
    '''
    Load All Reads
    '''
    print "Loading Reads..."
    SVM_BRF=sys.argv[2]
    BRF_lines=open(SVM_BRF).readlines()

    
    '''
    Locus-Based Filter
    '''
    print "Blacklist Filtering..."
    blacklist=build_blacklist(BRF_lines)
    BL_filtered_BRF_lines=run_blacklist(blacklist,BRF_lines)


    '''
    Read-Based Filter
    '''
    print "SVM Filtering..."
    THRESHOLD=0
    SVM_BL_filtered_BRF_lines=SVM_filter(BL_filtered_BRF_lines,THRESHOLD)
   

    '''
    R1R2 Filter
    '''
    print "R1R2 Filtering..."
    R1R2_blacklist=build_R1R2_blacklist(SVM_BL_filtered_BRF_lines)
    R1R2_SVM_BL_filtered_BRF_lines=run_R1R2(SVM_BL_filtered_BRF_lines, R1R2_blacklist)
    

    '''
    MAKE DETECTIONS
    '''
    VCF=sys.argv[1]
    tumor_dict=generate_tumor_dict(VCF)
    detected_reads=make_detections(R1R2_SVM_BL_filtered_BRF_lines,tumor_dict)
    reads_checked=count_fragments(BL_filtered_BRF_lines)
    sites_checked=count_sites(SVM_BL_filtered_BRF_lines)
    sites_detected=count_sites(detected_reads)
    detection_rate=sites_detected/float(reads_checked)

    '''
    WRITE REPORT
    '''
    BAM_id=sys.argv[2].split("_VS_")[0]
    VCF_id=sys.argv[2].split("_VS_")[1][:-8]
    VCF_id=VCF_id.split("_")[0]
    
    output_file=open(sys.argv[3],"a")
    output_file.write("BAM,VCF,sites checked, reads checked, sites detected, detection rate \n")
    output_file.write(",".join(map(str,[BAM_id, VCF_id, sites_checked, reads_checked, sites_detected, detection_rate,"\n"])))
    output_file.close()
main()
