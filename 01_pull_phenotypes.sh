# this script gets you phenotypes of you choice

################################################
########### COGNITIVE FUNCTION ###############
###############################################
#Category 100026
#NOTE: selection of best processable fields read out
############################################

# Response time on correctly identifying matches 20023
# Memory: max digits remebered 4282
# Fluid intelligence: score 20016
# Matrix pattern recognition
## Number of puzzles solved 6373
## Item slected for each puzzle 6332
## Duration spent answering the puzzle 6333
# Tower re-arranging: nb_puzzles_correct 21004
# Picture vocabulary: vocab_lvl 6364
# Symbol digit substitution: nb_sym_digit_match_corr 23324
# Word association: nb_corr_word_assocatio 20197
# Prospective memory: prosp_memory_result 20018
# Pair matching: nb_incorr_pair_matches_per_round 399

##########################################
## MEDICAL CONDICTIONS (Category 100044)###
##########################################
# Diastolic blood pressure, automated reading (Data-Field 4079)
## Two measures of blood pressure were taken a few moments apart.

# overall health rating (self-rated) 2178
## Excellent
## Good
## Fair
## Poor
## Do not know
## Prefer not to say

# Vascular/heart problems diagnosed by doctor 6150
## Heart Attack
## Angina
## Stroke
## High_BP
## None of the above
## Prefer not to say


#####################################
####### MEDICATIONS #################
#####################################
# Number of medications taken 137


# here is the code:
module load Python/3.13.1-GCCcore-14.2.0
module load SciPy-bundle/2022.05-foss-2022a
#source /cluster/projects/p33/users/maxk/py3/bin/activate

#python3 /cluster/projects/p33/users/maxk/UKB/ukb_helper.py pheno --input "/tsd/p33/data/durable/s3-api/ukblake/phenotypes/ukb*csv" --fields 31 33 34 52 53 54 --out /cluster/projects/p33/users/maxk/UKB/data/demo.txt

# 31 Sex, 33 date of birth, 34 year of birth, 52 month of birth
# 53	Date of attending assessment centre
# 54 UK Biobank assessment centre
# 129	Place of birth in UK - north co-ordinate	Early life factors
# 130	Place of birth in UK - east co-ordinate	Early life factors
# 1647 place of birth non-UK
# 20115 same as 1647?
# 738 Average total household income before tax
python3 /cluster/projects/p33/users/maxk/UKB/ukb_helper.py pheno --input "/tsd/p33/data/durable/s3-api/ukblake/phenotypes/ukb*csv" --fields 31 33 34 52 53 54 129 1647 130 129 20115 738 --out /cluster/projects/p33/users/maxk/UKB/environment/data/demo/environment


