process COVERAGE {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}/08_coverage/depth", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai), path(consensus)

    output:
    path "*.depth", emit: depth
    path "*.coverage.tsv", emit: stats
    path "versions.yml", emit: versions

    script:
    def prefix = meta.id
    """
    samtools depth -a -@ ${task.cpus} ${bam} > ${prefix}.depth
    total=\$(grep -v "^>" ${consensus} | tr -d '\\n' | wc -c)
    n_count=\$(grep -v "^>" ${consensus} | tr -d '\\n' | tr 'a-z' 'A-Z' | tr -cd 'N' | wc -c)
    percent_n=\$(awk -v n="\$n_count" -v len="\$total" 'BEGIN{printf "%.2f", (len>0?(n/len)*100:0)}')

    awk -v glen="\$total" -v sample="${prefix}" -v ncount="\$n_count" -v pctn="\$percent_n" '
    BEGIN{min=1e9;max=0;sum=0;sumsq=0;count=0;c1=c10=c20=c30=c50=0}
    {d=\$3;sum+=d;sumsq+=d*d;count++;if(d<min)min=d;if(d>max)max=d
     if(d>=1)c1++;if(d>=10)c10++;if(d>=20)c20++;if(d>=30)c30++;if(d>=50)c50++}
    END{
      mean=(count>0)?sum/count:0;std=(count>1)?sqrt((sumsq-(sum*sum)/count)/(count-1)):0
      if(min==1e9)min=0
      printf "Sample\\tMeanDepth\\tMinDepth\\tMaxDepth\\tStdDev\\tCov1x\\tCov10x\\tCov20x\\tCov30x\\tCov50x\\tN_count\\tGenome_length\\tPercent_N\\n"
      printf "%s\\t%.2f\\t%d\\t%d\\t%.2f\\t%.2f\\t%.2f\\t%.2f\\t%.2f\\t%.2f\\t%d\\t%d\\t%s\\n",
        sample,mean,min,max,std,(glen>0?(c1/glen)*100:0),(glen>0?(c10/glen)*100:0),
        (glen>0?(c20/glen)*100:0),(glen>0?(c30/glen)*100:0),(glen>0?(c50/glen)*100:0),ncount,glen,pctn
    }' ${prefix}.depth > ${prefix}.coverage.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
