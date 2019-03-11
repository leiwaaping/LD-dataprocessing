# LD-dataprocessing
genome data processing before LD analysis

 ## usage ##
**input** : target SNP file and vcf files(saved in one folder)
```
eg. snp.txt
rs2481074
rs75452934
rs6427562
rs6672793
rs7543952
rs2256505
rs2481073
rs12131275
rs7524408
rs2779800
```



**command** ：$ python LDprocessing.py -t <targetsnp> -i <vcfs_dir> -o <outputfile_name> -d <max_distance_between loci>
```
li@biointelligence:~$ python LDprocessing.py -t 02sampsnp.list -i vcf -o LD.txt -d 2000
target snp list in : 02sampsnp.list
vcfs saved in dir : vcf
output file : LD.txt
max distance between loci is : 2000
################### processing compro #####################
################### counting candidate loci #####################
Run time is 20 s

```



**output** ：output file marks info from two loci
 ```
 eg.LD.txt
SNP1 SNP2 #AB #A #a #B #b #Sample
rs9427129 rs9427406 2 2 4 3 3 5
rs9427129 rs945706 2 2 4 3 3 5
rs9427129 rs945707 1 2 4 3 3 5
rs9427129 rs945708 1 2 4 2 4 5
rs9427129 rs945710 2 2 4 3 3 5
rs9427129 rs945711 1 2 4 3 3 5
rs9427129 rs945712 2 2 4 2 4 5
rs9427129 rs945792 2 2 4 3 3 5
rs9427129 rs949836 1 2 4 2 4 5
```

also "tempunit.txt" file will be created at the same time, which records SNP genotype for each vcf.
```
rs      chr     position        ref             vcf/5.vcf       vcf/6.vcf       vcf/test.vcf    vcf/2.vcf       vcf/4.vcf       main_allele
rs9427123       1       165791783       G       A       G       G       G       A       A
rs9427127       1       165971604       C       T       T       C       C       T       T
rs9427128       1       165993765       C       C       T       C       C       T       T
rs9427129       1       166000228       G       A       G       G       A       A       A
rs9427406       1       162708770       G       A       G       G       G       A       A
rs945706        1       162076552       C       A       C       C       C       A       A
rs945707        1       162076770       T       T       C       T       T       C       C
rs945708        1       162118901       C       C       T       C       T       T       T
rs945710        1       162139270       T       C       T       T       T       C       C

```
