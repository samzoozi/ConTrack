# ConTrack
ConTrack is a scalable concept tracking model, which can track multiple concepts and participation of multiple actors in each concept over time, in high dimensional temporal data. 

The detailed information about the model is provided in the following paper:
Zonoozi, Ali, et al. "ConTrack: A Scalable Method for Tracking Multiple Concepts in Large Scale Multidimensional Data." Data Mining (ICDM), 2016 IEEE 16th International Conference on. IEEE, 2016.

 
 ## Framework
The model is implemented in Java and runs on "Petuum" distributed framerwork (https://www.petuum.com/), which requires Gradle (https://gradle.org) to be installed as the automation and compilation toolkit.

For more information please refer to JBösen library (https://github.com/sailing-pmls/jbosen). JBösen is the Java implementation of the Petuum project which allows applications to run in parallel using a parameter server architecture. The contents of this project should be placed under a new directory under the 'app' directory, and the name of the parent directory has to be included in the settings.gradle file.
