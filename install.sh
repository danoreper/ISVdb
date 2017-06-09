echo "Copying files from Zenodo"
curl https://zenodo.org/record/399474/files/File_S2.tar.gz > ./data/cc_probs_2016-03_mm10.tar.gz
cd ./data/
tar -xzvf cc_probs_2016-03_mm10.tar.gz
cd -

curl  https://zenodo.org/record/399474/files/File_S3.gz > ./data/vcf/rel1410/mgp.v4.snps.dbSNP.vcf.gz

curl  https://zenodo.org/record/399474/files/File_S4.gz > ./data/vcf/rel1410/mgp.v4.indels.dbSNP.vcf.gz

curl  https://zenodo.org/record/399474/files/File_S5.gz > ./data/vcf/rel1410/mgp.v4.snps.MT.dbSNP.vcf.gz


curl https://zenodo.org/record/399474/files/File_S7.gtf.tar.gz > ./data/b6_reference/mus_musculus.GRCm38.75/Mus_musculus.GRCm38.75.gtf.tar.gz
cd ./data/b6_reference/mus_musculus.GRCm38.75
tar -xzvf Mus_musculus.GRCm38.75.gtf.tar.gz  -O > Mus_musculus.GRCm38.75.gtf
cd -

