version 1.0

workflow xin_shai_panel {
    
    input {
        String inputSamplesFile
        Array[Array[File]] inputSamples = read_tsv(inputSamplesFile)

        String ref
        String outdir

        String java_path
        String picard_path
        String samtools_path 
        String bwa_path
        String varscan_path
        String snpEff_path
    }

  scatter (sample in inputSamples) {
    call run_fastp {
      input:  
        sample_id=sample[0],
        outdir= outdir + "/untrim/fastp/" + sample[0] + "/", 
        fastq1=sample[1],
        fastq2=sample[2],
    }
    call run_align as run_align1  {
        input:
            fastq1=run_fastp.fq_out1,
            fastq2=run_fastp.fq_out2,
			sample_id=sample[0],
            outdir=outdir + "/untrim/sortSam/" + sample[0] + '/',
            ref=ref,
            bwa_path=bwa_path,
    }

    call run_SortSam as run_SortSam1  {
        input:
            java_path=java_path,
            picard_path=picard_path,
            outdir=outdir + "/untrim/sortSam/" + sample[0] + '/',
        sample_id=sample[0],
            sam_file=run_align1.align_sam_out,
            samtools_path = samtools_path, 
    }

    call trimPrimer2fq  {
        input:
        sample_id=sample[0],
            in_bam=run_SortSam1.bam_out,
            in_bam_bai=run_SortSam1.bam_out_bai,
            outdir=outdir + '/untrim/trimPrimer/' + sample[0] + '/',
    }

# 4 重复对新数据再跑一次

    call run_align  as run_align2  {
        input:
            fastq1=trimPrimer2fq.out_fq1,
            fastq2=trimPrimer2fq.out_fq2,
        sample_id=sample[0],
            outdir=outdir + "/clean_trim/sortSam/"  + sample[0] + '/',
            ref=ref,
            bwa_path=bwa_path,
    }

    call run_SortSam as run_SortSam2  {
        input:
            java_path=java_path,
            picard_path=picard_path,
            outdir=outdir + '/clean_trim/sortSam/' + sample[0] + '/',
        sample_id=sample[0],
            sam_file=run_align2.align_sam_out,
            samtools_path = samtools_path, 
    }


    call run_varscan {
        input :
            java_path=java_path,
            outdir=outdir + '/clean_trim/varscan/'+ sample[0] + '/',
            in_sort_bam = run_SortSam2.bam_out,  
            in_sort_bam_bai = run_SortSam2.bam_out_bai,  
        	sample_id=sample[0],
            samtools_path = samtools_path, 
            snpEff_path =snpEff_path ,
            varscan_path =varscan_path,
            ref=ref,
     }

#    call reohit2   {
#        input:
#            outdir=outdir + '/clean_trim/reohit2/'+ sample[0],
#        	sample_id=sample[0],
#            in_vcf=run_varscan.out_vcf,
#    }
  }

}


task run_fastp {
    input {
        String sample_id
        String fastp_path
        # reads pair
        File fastq1
        File fastq2
        String outdir
        # report parameters
        String html_report 
        String json_report 
    }
    String fastp_dir = outdir 
    String fq_out1 = fastp_dir +  sample_id + "_1.clean.fq.gz"
    String fq_out2 = fastp_dir +  sample_id + "_2.clean.fq.gz"
    command <<<
        mkdir  -p ~{fastp_dir} 
        ~{fastp_path}  \
        -i ~{fastq1} \
        -I ~{fastq2} \
        -o ~{fq_out1} \
        -O ~{fq_out2} \
        -h ~{fastp_dir}/~{sample_id}.fastp.html \
        -j ~{fastp_dir}/~{sample_id}.fastp.json 

    >>>


output {
        File fq_out1 = fq_out1
        File fq_out2 = fq_out2
    }
}



task run_align {
    input {
        File fastq1
        File fastq2
        String outdir
        String sample_id
        String ref

        String opts = " -t 6 -M "
                             #"'@RG\\tID:FH_B_30\\tPL:illumina\\tPU:FH_B_30\\tLB:FH_B_30\\tSM:FH_B_30\\tCN:REO'"
        String reads_group =  "'@RG\\tID:" + sample_id + "\\tPL:illumina\\tPU:"+ sample_id +"\\tLB:" + sample_id +"\\tSM:"+ sample_id +"\\tCN:REO'"

        String bwa_path 


        # Resource
        String memory = "32G"
        
    }

    String sam_dir = outdir  
    String sam_out = sam_dir + sample_id + ".sam"

    command <<<
        mkdir -p ~{sam_dir}
        ~{bwa_path} mem ~{opts} \
        -R ~{reads_group} \
         ~{ref} \
         ~{fastq1} \
         ~{fastq2} \
        > ~{sam_out} \
    >>>

    runtime {
        memory: memory
    }

    output {
        File align_sam_out = sam_out
    }
}

task run_SortSam {
    input {
        String outdir
        String sample_id

        String java_path 
        String picard_path 
        String sam_file 
        String samtools_path 
        String opts = "  -Xmx3g   -XX:MaxPermSize=512m -XX:-UseGCOverheadLimit "

        # Resource
        String memory = "32G"
        
    }

    String bam_out = outdir  + sample_id + ".sort.bam"
    String bam_out_bai = outdir + "/" + sample_id + ".sort.bam.bai"

    command <<<
        mkdir -p ~{outdir}
        ~{java_path}  ~{opts} \
         -jar ~{picard_path} \
         SortSam \
         I=~{sam_file} \
         O=~{bam_out} \
         SO=coordinate

         ~{samtools_path} index ~{bam_out}
        /bin/rm   ~{sam_file}
    >>>

    runtime {
        memory: memory
    }

    output {
        File bam_out = bam_out
        File bam_out_bai = bam_out_bai
    }
}


task trimPrimer2fq {
    input {
        String sample_id
        File in_bam
        String bam_trimPrimer2fq_path

        String outdir
        File in_bam
        File in_bam_bai
        # Resource
        String memory = "32G"
        
    }

    String fq_path = outdir 
    String fq_out1 = fq_path+"/" +sample_id + "_1.fq.gz"
    String fq_out2 = fq_path+"/" +sample_id + "_2.fq.gz"

    command <<<
        mkdir  -p ~{fq_path}
        python  ~{bam_trimPrimer2fq_path} \
         ~{in_bam} \
         ~{fq_path}/~{sample_id} 
        gzip -f   ~{fq_path}/~{sample_id}_1.fq
        gzip -f  ~{fq_path}/~{sample_id}_2.fq
    >>>

    runtime {
        memory: memory
    }

    output {
        File out_fq1 = fq_out1
        File out_fq2 = fq_out2
    }
}



task run_varscan {
    input {
        String outdir
        String sample_id

        String in_sort_bam 
        String in_sort_bam_bai 
        String java_path 
        String samtools_path 
        String snpEff_path 
        String varscan_path
        String ref
    }

    String varscan_results_path = outdir  
    String out_vcf = varscan_results_path + "/" +sample_id+".varscan.snpIndel.anno.vcf.gz"  

    command <<<

        mkdir -p  ~{varscan_results_path}

        ~{samtools_path} mpileup -f ~{ref} ~{in_sort_bam}  | ~{varscan_path} mpileup2cns --min-coverage 4   --output-vcf 1 > ~{varscan_results_path}/~{sample_id}.varscan.cns.vcf &&\
            egrep -v '0/0|\./\.' ~{varscan_results_path}/~{sample_id}.varscan.cns.vcf > ~{varscan_results_path}/~{sample_id}.varscan.snpIndel.vcf && \
            gzip -f  ~{varscan_results_path}/~{sample_id}.varscan.cns.vcf
            gzip -f ~{varscan_results_path}/~{sample_id}.varscan.snpIndel.vcf && \
            ~{java_path}  -jar ~{snpEff_path} -noShiftHgvs -nextProt -lof -oicr GRCh38.p14.RefSeq.2209 ~{varscan_results_path}/~{sample_id}.varscan.snpIndel.vcf.gz > ~{varscan_results_path}/~{sample_id}.varscan.snpIndel.anno.vcf && \
            gzip -f  ~{varscan_results_path}/~{sample_id}.varscan.snpIndel.anno.vcf 
    >>>

#    runtime {
#        memory: memory
#    }

    output {
        File  out_vcf = out_vcf 
    }
}


task reohit2 {
    input {
        String python2_path
        String reohit2_path
        String reohit2_db
        String xinshai_bed
        String in_vcf
        String sample_id

        String outdir
        # Resource
        String memory = "32G"
        
        String out_path = outdir + "/" + sample_id +'.cvf'
    }


    command <<<
    mkdir -p ~{outdir}
    ~{python2_path} ~{reohit2_path} --database ~{reohit2_db}  -V ~{in_vcf}  -N 1 -B 13 -O ~{out_path}  -E ~{xinshai_bed} 

    >>>

    runtime {
        memory: memory
    }

    output {
#File out_fq2 = fq_out2
    }
}

