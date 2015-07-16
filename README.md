# Pre_ICA_Cleaning
Uses median windowed mutual information reduction index to find the best method for cleaning EEG data before performing ICA.

The main function is pre_ica_cleaning_index(). It calculates an index (Pre ICA Cleaning Index, or PICI) that can be used to compare the quality of ICA cleaning (e.g. artifact rejection) using median windowed mutual infromation method. The higher this value is, the better is the pre-ica cleaning method.

In this method, Mutual information reduction  index, or MIR is calculated in short ~1 s windows of the original (full, non rejected or otherwise cleaned data) and the median of these MIR values is calculated for both an ICA decomposition (often calculated on a separate, clean data) and an sphering matrix calculated based on a geometric-median calculated covariance matrix.

To use this function you should run the same ICA algorthms (e.g Infomax) on the 'clean' portion of the data, and then use provide this ICA matrix along with FULL NON-REJECTED DATA in EEG input variable to this function. DO NOT send cleaned data in nonCleanedEEG variable to this function as this would certainly produce misleading results. The cleaning method that produces the highest PICI value is more suitable to be used in later analysis. Bad channels should be removed from nonCleanedEEG and all other EEG datasets before using this function.
