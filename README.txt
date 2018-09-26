1. Put the four files
"collector.R", "tester.R", "trainingset.csv", "trees.dat"
in your working directory.

2. From your R console, type
source("tester.R")
to run the testing program.

3. Select the .csv file you want to test.

4. If trainingset.csv was updated prior to the latest run, enter "Y" when prompted "Update Model?".
Using the already existing model reduces running time.

5. The prediction details is wrote to "predictiondata.csv". Visualization of the aggregate prediction is wrote to "timeline.png".

6. Run "collector.R" to add a training data to the existing training set.