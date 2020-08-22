# RMO-based algorithm
The RMO-based algorithm could be used to identify the origin of brain metastasis tumors by the within-sample relative methylation orderings (RMOs) of the CpG sites.
# Usage
selec.stable.pair(Beta,cut_off)

Input
|Arguments|Description|
|---------|-----------|
|Beta|The DNA methylation profiles of samples of a specific cancer type. The row name of methylation profiles was CpG probe name, and the column name of methylation profiles was sample ID (or sample name).|
|cut_off|The criteria for identifying stable CpG site pairs in this cancer type. The default setting of freq is 0.95.|

Output
|Description|
|-----------|
|The output file was stable CpG site pairs for the analyzed cancer type. The output file was a data.frame with two columns. For each row, the methylaton level of the CpG site in the first column was stably higher than that of the GpG site in the second column.|

classification(reversal_pair,Beta,reversal_pair_cancer)
Input
|Arguments|Description|
|---------|-----------|
|reversal_pair|The reversal CpG site pairs. A list of reversal CpG site pairs was a data.frame. In this study, seven lists of reversal CpG site pairs between GM and each of the other cancer types will be integrated into a list format. And each data.frame (or element) in this list needs to be named this form - 'cancer-cancer pair', such as GM-BRC.|
|Beta|The DNA methylation profiles of testing samples. The row name of methylation profiles was CpG probe name, and the column name of methylation profiles was sample ID (or sample name).|
|reversal_pair_cancer|A cancer type. For the input reversal_pair, the methylaton level of the CpG sites in the first column are stably higher than that of the corresponding GpG sites in the second column, which occurred in this cancer type.| 

Output
A list of three elements
|Name|Description|
|---------|-----------|
|Output$prediction_matrix|A matrix of classification result of each sample in each reversal pair. If only one inversion pair is input, the result is the same as Output$prediction_label.|
|Output$prediction_label|Classification result of each sample.|
|Output$existence_reversal_pair|A list of reversal pairs. It outputs the reversal pairs available in the input methylation profiles.|

calculate.similarity(reversal_pair,Beta,reversal_pair_cancer)
Input
|Arguments|Description|
|---------|-----------|
|reversal_pair|The reversal CpG site pairs. A list of reversal CpG site pairs was a data.frame. In this study, seven lists of reversal CpG site pairs between GM and each of the other cancer types will be integrated into a list format. And each data.frame (or element) in this list needs to be named this form - 'cancer-cancer pair', such as GM-BRC.|
|Beta|The DNA methylation profiles of testing samples. The row name of methylation profiles was CpG probe name, and the column name of methylation profiles was sample ID (or sample name).|
|reversal_pair_cancer|A cancer type. For the input reversal_pair, the methylaton level of the CpG sites in the first column are stably higher than that of the corresponding GpG sites in the second column, which occurred in this cancer type.| 

Output
A list of three elements 
|Name|Description|
|---------|-----------|
|Output$similarity|A matrix of similarity scores to tumor of each sample.|
|Output$prediction_label|The cancer type with the highest similarity score for each testing sample.|
|Output$existence_reversal_pair|A list of reversal pairs. It outputs the reversal pairs available in the input methylation profiles.|
