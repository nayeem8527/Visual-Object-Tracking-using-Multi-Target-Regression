This is a preliminary implementation of the Online Discriminative Dictionary Learning (ODDL) tracking method proposed in 

"Online Discriminative Dictionary Learning for Visual Tracking",
Fan Yang, Zhuolin Jiang and Larry S. Davis,
IEEE Winter Conference on Applications of Computer Vision (WACV), pp. 854-861, 2014.

It has been tested with MATLAB R2012a on Windows 7/8 and MATLAB R2013b on Ubuntu 14.04. It should also work well on other platforms, provided that you have successfully compiled required external libraries.

HOW TO RUN

1. Specify the parameters in "trackparam". You should keep all parameters fixed except the sampling parameters for individual sequences.

2. Change the directory of test sequences in "tracker": seq_path. A sample sequences is included in "data" directory.

3. run "tracker" to see how the tracker proceeds. Results (figures and locations of bounding boxes) are saved in "results" directory.

DEPENDENCIES

This implementation requires a few external libraries.
- VLfeat to extract HoG features (http://www.vlfeat.org/)
- KSVDBox and OMPBox (http://www.cs.technion.ac.il/~ronrubin/software.html)
- Sparse Modeling Software (http://spams-devel.gforge.inria.fr/)

CONTACT

Please contact Fan Yang (fyang@umiacs.umd.edu) for any questions.