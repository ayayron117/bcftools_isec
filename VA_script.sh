#!/bin/bash

# Two ancestral samples: 3KO_611 & 3KO_674
# Four parental samples: 3KO_850 & 3KO_853, and 3KO_864 & 3KO_865
# Two surviving offspring: 3KO_894 and 3KO_900

module load bcftools/1.17

# Path to annotated vcf files in sarek results folder
annots=$(pwd)/sarek_output/annotation/haplotypecaller

# Create new directory that will contain results from this analysis.
# It will be called "variant_analysis"
main=$(pwd)/variant_analysis
mkdir -p $main

# Create a sub directory in the variant_analysis folder that will contain
# copies of the annotated vcf files for the ancestors. Copy those files to
# that directory.
ancestors=$main/ancestors
mkdir -p $ancestors
cp -r $annots/3KO_611_S14/3KO_611_S14* $ancestors
cp -r $annots/3KO_674_S2/3KO_674_S2* $ancestors

# Create a sub directory in the variant_analysis folder called "group_A" and then
# create another sub directory in this folder called "parents" that will contain 
# the annotated vcf files of the parents of one of the offspring. Also create 
# another sub directory in group_A called "offspring" that contains the file for the 
# offspring of those parents. Copy the files to their respective directories.
GA=$main/group_A
parents_A=$GA/parents
offspring_A=$GA/offspring
mkdir -p $GA
mkdir -p $parents_A
mkdir -p $offspring_A
cp -r $annots/3KO_850_S3/3KO_850_S3* $parents_A
cp -r $annots/3KO_853_S4/3KO_853_S4* $parents_A
cp -r $annots/3KO_894_S1/3KO_894_S1* $offspring_A

# Do the same for the other group
GB=$main/group_B
parents_B=$GB/parents
offspring_B=$GB/offspring
mkdir -p $GB
mkdir -p $parents_B
mkdir -p $offspring_B
cp -r $annots/3KO_864_S5/3KO_864_S5* $parents_B
cp -r $annots/3KO_865_S6/3KO_865_S6* $parents_B
cp -r $annots/3KO_900_S2/3KO_900_S2* $offspring_B

# Merge the ancestor vcf files into one.
bcftools merge -m all -Oz -o $ancestors/ancestors_merged.vcf.gz \
${ancestors}/3KO_611_S14.haplotypecaller.filtered_snpEff_VEP.ann.vcf.gz \
${ancestors}/3KO_674_S2.haplotypecaller.filtered_snpEff_VEP.ann.vcf.gz

# Index the merged vcf file
bcftools index --tbi $ancestors/ancestors_merged.vcf.gz

####################################################################
###### Offspring variants in both parents & not in ancestors #######
####################################################################

# Use a for loop to get the intersection of variants in both parents of each group.
for parents in $parents_A $parents_B; do
    
    # Change working directory to the parents folder
    cd $parents
    
    # Create a vector that contains the file names of the parent files. We are using the
    # list command to do this.
    parent_files=($(ls *.vcf.gz))
    # echo "${parent_files[0]}"
    # echo "${parent_files[1]}"
    # Extract IDs from filenames to create output directory and file name. We are doing 
    # this by using a combination of the echo and cut commands. We need to use the echo
    # command because the cut command only takes in files as input. The output from using
    # echo will be treated as a file when piped into the cut command. The -d option in cut
    # stands for deliminator. The -f option stands for fields. We are setting the deliminator
    # to be "_", which means the file names contained in the parent_files vector will be split
    # based on the underscore symbol. Think of it as replacing the underscores in the file
    # name with spaces. The words left after removing the underscores are the fields.
    # Here's an example for one of the files: 
    # 3KO 850 S3.haplotypecaller.filtered snpEff VEP.ann.vcf.gz
    # I set this up so that the second field will be extracted from the file name.
    # In the example above, "850" is what will be extracted from that file name. 
    # This is done for each of the parent files and each ID that is extracted is stored in 
    # variables called ID1 and ID2.
    ID1=$(echo "${parent_files[0]}" | cut -d '_' -f 2)
    ID2=$(echo "${parent_files[1]}" | cut -d '_' -f 2)
    # echo "$ID1"
    # echo "$ID2"

    # Create a variable called fname that will be used to set the file name of the vcf
    # file that will contain the intersected variants of the parents. It will be a 
    # combination of ID1 and ID2.
    fname="${ID1}_${ID2}_intersection.vcf.gz"
    # echo "$fname"

    # Also use those IDs to create variable called "dir_name" that will contain the 
    # name of the directory that will be designated for the output of using bcftools isec. 
    dir_name="${ID1}_${ID2}_isec"
    # echo "$dir_name"

    # Perform the intersection using bcftools isec and use the -p to set the name of 
    # the output directory, which in this case is the name contained in "dir_name".
    # For the input files, use the 2 file names contained in the "parent_files" vector.
    # "-Oz" means save the vcf files in .gz format. 
    bcftools isec -p $dir_name -Oz ${parent_files[0]} ${parent_files[1]}
    
    # Copy the vcf file in the output folder that contains the intersection to the current working directory 
    # (remember, we are currently in one of the parent folders). Bcftools names this file "0002.vcf.gz".
    # Below are what each of the output files are:
    # 0000.vcf: Variants unique to first parent
    # 0001.vcf: Variants unique to second
    # 0002.vcf: Variants common to both files (the intersection) 
    # After copying the file, we also index it in .tbi format.
    cp $dir_name/0002.vcf.gz ./$fname
    bcftools index --tbi $fname

done

# Change working directory to the "variant_analysis" folder.
cd $main

# Get the intersection of variants in both parents and offspring using the parent intersection files
# and their respective offspring of each group.
for group in $GA $GB; do
    
    # Change working directory to one of the groups 
    cd $group

    # Change working directory to the parent folder of that group and get the name of the parent 
    # intersection file.
    cd parents
    parents_int_file=$(ls *_intersection.vcf.gz)

    # Change working directory to the offspring folder of the group and get the name of the offspring file.
    cd ../offspring
    offspring_file=$(ls *.vcf.gz)
    
    # Extract the IDs from the variables containing the filenames and use them to create the name of 
    # the output directory that will contain the output from bcftools isec. Also use them to set the name
    # of the intersection file.
    ID1=$(echo ${parents_int_file} | cut -d '_' -f 1)
    ID2=$(echo ${parents_int_file} | cut -d '_' -f 2)
    ID3=$(echo ${offspring_file} | cut -d '_' -f 2)
    fname="${ID1}_${ID2}_${ID3}_intersection.vcf.gz"
    dir_name="${ID1}_${ID2}_${ID3}_isec"

    # Change the working directory back to the group directory and create a folder called "homo" there.
    cd $group
    mkdir -p homo
    

    # Perform the intersection using bcftools isec. Save the output in the homo directory within a folder
    # that will be named using "dir_name". Use the paths to the parent intersection file and the offspring
    # file as input.
    bcftools isec -p homo/$dir_name -Oz parents/$parents_int_file offspring/$offspring_file
    
    # Copy the the output file that contains the intersection to the homo directory and index it.
    cp ./homo/$dir_name/0002.vcf.gz ./homo/$fname
    bcftools index --tbi ./homo/$fname

    # Use bcftools isec to determine the unique variants that are in both parents and offspring, but
    # are not in the ancestors. The file that will contain that will be the 0000.vcf.gz file.
    # Copy it to the "homo" directory and call it "homo_variants_not_in_ancestors.vcf.gz"
    bcftools isec -p homo/descendents_ancestors_isec -Oz ./homo/$fname $ancestors/ancestors_merged.vcf.gz
    cp homo/descendents_ancestors_isec/0000.vcf.gz ./homo/homo_variants_not_in_ancestors.vcf.gz
    bcftools index --tbi ./homo/homo_variants_not_in_ancestors.vcf.gz

done

#################################################################################
####### Offspring variants in either parent (not both) & not in ancestors #######
#################################################################################

# Change working directory to the "variant_analysis" folder.
cd $main

# For each group do the following:
# 1. Get the intersection between one parent and the offspring
# 2. Get the intersection between the other parent and the offspring
# 3. Run bcftools isec using those two files
# 4. Merge the vcfs that contain the unique variants from each file
# 5. Run isec using that merged file and the merged ancestors file
# 6. Re-name the resulting 0000.vcf.gz file 
# 
# Example for group A:
# 1. Intersection between 3KO_850 & 3KO_894 -> 850_894_intersection
# 2. Intersection between 3KO_853 & 3KO_894 -> 853_894_intersection
# 3. isec 850_894_intersection 853_894_intersection
# 4. merge 0000.vcf 0001.vcf
# 5. isec merged.vcf ancestors.vcf
# 6. Re-name 0000.vcf

for group in $GA $GB; do
    
    # Save the name of the group in a variable
    # ($group contains the full path to the group folder that the for loop 
    # is acting on, so we're using the command "basename" to extract the string 
    # that follows the last "/" in the full path)
    gname=$(basename "$group")

    # Change working directory to the parents folder
    cd $group/parents
    
    # Create a vector that contains the file names of the parent files.
    parent_files=($(ls *.vcf.gz))
    # echo "${parent_files[0]}"
    # echo "${parent_files[1]}"

    # Extract the IDs for each parent.
    ID1=$(echo "${parent_files[0]}" | cut -d '_' -f 2)
    ID2=$(echo "${parent_files[1]}" | cut -d '_' -f 2)
    # echo "$ID1"
    # echo "$ID2"

    # Change working directory to offspring folder
    cd $group/offspring

    # Extract the ID of the offspring
    offspring_file=$(ls *.vcf.gz)
    ID3=$(echo ${offspring_file} | cut -d '_' -f 2)

    # Create variables that will contain the names for each of the output directories and files
    dir_name1="${ID1}_${ID3}_isec"
    fname1="${ID1}_${ID3}_intersection.vcf.gz"

    dir_name2="${ID2}_${ID3}_isec"
    fname2="${ID2}_${ID3}_intersection.vcf.gz"
 
    # Change working directory to the group directory
    cd $group

    # Create a folder called hetero and cd into it
    mkdir -p hetero
    cd hetero

    # Perform the intersections using bcftools isec
    bcftools isec -p $dir_name1 -Oz ../parents/${parent_files[0]} ../offspring/$offspring_file
    bcftools isec -p $dir_name2 -Oz ../parents/${parent_files[1]} ../offspring/$offspring_file
    
    # Copy the intersection files to the hetero directory and re-name them
    cp $dir_name1/0002.vcf.gz ./$fname1
    cp $dir_name2/0002.vcf.gz ./$fname2

    # Index the files
    bcftools index --tbi $fname1
    bcftools index --tbi $fname2

    # Run isec using both of the intersection files
    dir_name3="${ID1}_${ID3}-${ID2}_${ID3}_isec"
    bcftools isec -p $dir_name3 -Oz $fname1 $fname2

    # Merge the vcf files that contain the unique variants from each file
    # (Merge 0000.vcf.gz and 0001.vcf.gz)
    bcftools merge -m all -Oz -o ${gname}_hetero.vcf.gz \
    $dir_name3/0000.vcf.gz \
    $dir_name3/0001.vcf.gz

    # Index
    bcftools index --tbi ${gname}_hetero.vcf.gz

    bcftools isec -p descendents_ancestors_isec -Oz ${gname}_hetero.vcf.gz $ancestors/ancestors_merged.vcf.gz
    cp descendents_ancestors_isec/0000.vcf.gz ${gname}_hetero_vars_not_in_ancestors.vcf.gz
    bcftools index --tbi ${gname}_hetero_vars_not_in_ancestors.vcf.gz

done