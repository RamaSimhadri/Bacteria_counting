# Bacteria_counting
Counting the number of bacteria in a confocal Z-stack

This script will count the number of Wolbachia bacteria in a confocal Z-stack. 
FIJI or ImageJ have scripts for edge detection and counting objects, but this script does better when the bacteria are clumped together. 
Once the clumps are identified, local maxima, i.e. the centers of each bacterium within the clump are counted, which produces a more accurate estimate of the number of bacteria in the image.

Input: Confocal Z-stacks of cell lines. 
copy the script into the directory with the TIF images of the cofocal sections of your bacteria or embed this in a wrapper to go to the directory of interest

Output: Plots of the bacteria identified and the number of bacteria in each plane. 
The list of bacteria can be exported into an excel file or CSV file using xlswrite or csvwrite
