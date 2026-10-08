# Step 9 Biocorrelation Details

In the `metadata` table the columns with a bioactivity score must be named with the prefix "BIOACTIVITY_" and the columns 
with the correlation groups (samples to be used in each correlation) must be named with the prefix "COR_". All pairwise
combinations of bioactivity scores and correlation groups will create a new biocorrelation score store in a different 
column of the count table named as "COR_<BIOACTIVITY suffix>_<COR suffix>".

The correlation score will produce NA values with warnings if the counts of a spectra (quantification) in the selected 
samples are all equal 0 or if the bioactivity of the selected samples are all the same. 
And, it will produce "CTE" if the counts of the selected samples have the same values (standard deviation equals 0).