# User data aggregation for StallCatchers.com Crowdsourcing Initiative

## 1. INTRODUCTION 
The goal of crowdsourcing and machine learning is to combine the human knowledge and machine computing to aid research and problem solving in areas where combined efforts of human and machine are greatly needed. Crowdsourcing creates new training data for machine learning. For instance, Zhu et al.[1] used Amazon Mechanical Turk [2] to acquire ground truth for their generative adversarial networks. However, crowdsourcing faces a general problem that expertise of the majority could be questionable. When the crowd makes contradictory predictions or even sometimes crowd is wrong [3], it is a challenge for scientists to determine the optimal result. Therefore, it is necessary to develop an algorithm to adaptively combine crowd results and make a decision based on individual’s precisions from previously known training data. 
## 2. DATASET 
### 2.1 Data acquisition 
Data were acquired through an online game stallcatchers.com, which has more than 4000 users. Players start with a training video to learn how to play the game followed by a hands-on session that after each round of the game, they receive details feedback describing the correct answer. Each player at each round of the game sees a 3D stack image of brain blood vessels, which is acquired using multiphoton microscopy. The targeted vessel was outlined for users to determine whether the vessel was stalled or flowing. Due to rastering mechanism of multiphoton microscopy, different slices were acquired within half a second. This fact allows users to observe the movement of red blood cells as black patches and help them to discretize between stalled and flowing vessels. In order to record user accuracy, after several novel images, they will be tested on one of the calibration vessels with known answer. 
### 2.2 Dataset structure 
The StallCatchers data is managed using MySQL. The users data can be exported as csv text files that each line of record contain timestamp, vessel ID, user ID, user answer (0 = flowing vessel, 1 = stalled vessel), and Calibration flag (0 = novel vessel, 1 = calibration vessel). 

#### REFERENCE 
[1] J.-Y. Zhu, T. Park, P. Isola, and A. A. Efros, “Unpaired Image-to-Image Translation using Cycle-Consistent Adversarial Networks,” arXiv, 30-Mar-2017. 

[2] “Amazon Mechanical Turk - Welcome.” [Online]. Available: https://www.mturk.com. [Accessed: 19-Sep-2017]. 

[3] D. Prelec, H. Sebastian Seung, and J. McCoy, “A solution to the single-question crowd wisdom problem,” Nature, vol. 541, no. 7638, pp. 532–535, 2017. 

[4] A. Beygelzimer, S. Dasgupta, and J. Langford, “Importance weighted active learning,” in Proceedings of the 26th Annual International Conference on Machine Learning - ICML ’09, 2009. 

[5] M. Joglekar, H. Garcia-Molina, and A. Parameswaran, “Comprehensive and reliable crowd assessment algorithms,” in 2015 IEEE 31st International Conference on Data Engineering, 2015. 

[6] J. Whitehill, T.-F. Wu, J. Bergsma, J. R. Movellan, and P. L. Ruvolo, “Whose Vote Should Count More: Optimal Integration of Labels from Labelers of Unknown Expertise,” in Advances in Neural Information Processing Systems, 2009, pp. 2035–2043. 

[7] V. C. Raykar et al., “Supervised learning from multiple experts,” in Proceedings of the 26th Annual International Conference on Machine Learning - ICML ’09, 2009. 

[8] S. K. Warfield, K. H. Zou, and W. M. Wells, “Simultaneous truth and performance level estimation (STAPLE): an algorithm for the validation of image segmentation,” IEEE Trans. Med. Imaging, vol. 23, no. 7, pp. 903–921, Jul. 2004. 
