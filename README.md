# Visual-Object-Tracking-using-Multi-Target-Regression

Proposing a visual tracking technique using Multi-layer Multi-target regression(MMR) which enables simultaneous modelling of inter-target correlations and input-output relationship via robust low-rank learning algorithm. The proposed method takes in the fhog features from the current frame and with the help of regression we classify the region as foreground or background. The regression model is updated in an online fashion after every fifth frame to capture the change in frames for efficient tracking.
For more detail explanation please read [here](https://drive.google.com/file/d/1dW2MneJdd1d0lOsc67kJ5TUFbclIGE8I/view?usp=sharing).

# Results

**Tracking Qualitative Results**

|![ball1_global_0009](https://user-images.githubusercontent.com/10145585/50384238-53d29600-06e8-11e9-9d27-e0c367b4178f.png)|![ball1_global_0046](https://user-images.githubusercontent.com/10145585/50384244-6a78ed00-06e8-11e9-83c2-6432f00dafb3.png)|![ball1_global_0079](https://user-images.githubusercontent.com/10145585/50384266-9a27f500-06e8-11e9-9cc2-054886a03874.png)|![ball1_global_0097](https://user-images.githubusercontent.com/10145585/50384270-a318c680-06e8-11e9-953e-f7b088035a1f.png)|![ball1_global_0105](https://user-images.githubusercontent.com/10145585/50384274-ac099800-06e8-11e9-9baf-287ae5a8a677.png)
|:---:|:---:|:---:|:---:|:---:|

|![bag_global_0001](https://user-images.githubusercontent.com/10145585/50384364-eb84b400-06e9-11e9-9f6c-611068a9e415.png)|![bag_global_0037](https://user-images.githubusercontent.com/10145585/50384367-fa6b6680-06e9-11e9-826b-66341ddf32dc.png)|![bag_global_0075](https://user-images.githubusercontent.com/10145585/50384372-0bb47300-06ea-11e9-980a-49767138b2e8.png)|![bag_global_0137](https://user-images.githubusercontent.com/10145585/50384375-15d67180-06ea-11e9-84af-5359bda85d80.png)|![bag_global_0187](https://user-images.githubusercontent.com/10145585/50384381-1ec74300-06ea-11e9-9b21-0ed0247118bc.png)
|:---:|:---:|:---:|:---:|:---:|

|![tiger1_global_0009](https://user-images.githubusercontent.com/10145585/50384391-3dc5d500-06ea-11e9-924c-317740b24312.png)|![tiger1_global_0097](https://user-images.githubusercontent.com/10145585/50384397-49b19700-06ea-11e9-8e22-2399949ad424.png)|![tiger1_global_0198](https://user-images.githubusercontent.com/10145585/50384402-5209d200-06ea-11e9-8155-1751278ff995.png)|![tiger1_global_0314](https://user-images.githubusercontent.com/10145585/50384406-5a620d00-06ea-11e9-85c2-e3f7528395f6.png)|![tiger1_global_0352](https://user-images.githubusercontent.com/10145585/50384410-63eb7500-06ea-11e9-8764-1e0d8f7337de.png)
|:---:|:---:|:---:|:---:|:---:|

|![surfer_global_0009](https://user-images.githubusercontent.com/10145585/50384416-72d22780-06ea-11e9-9483-5216841f0dec.png)|![surfer_global_0041](https://user-images.githubusercontent.com/10145585/50384419-7d8cbc80-06ea-11e9-85d8-5a788a94728f.png)|![surfer_global_0142](https://user-images.githubusercontent.com/10145585/50384422-8aa9ab80-06ea-11e9-85a6-7fba43cbdaf6.png)|![surfer_global_0243](https://user-images.githubusercontent.com/10145585/50384426-95644080-06ea-11e9-878b-a4db9d3bfcd9.png)|![surfer_global_0338](https://user-images.githubusercontent.com/10145585/50384483-2d622a00-06eb-11e9-8292-c4501a873c12.png)
|:---:|:---:|:---:|:---:|:---:|

|![boy_global_0016](https://user-images.githubusercontent.com/10145585/50384448-b9278680-06ea-11e9-80eb-e7415aa3ad5c.png)|![boy_global_0055](https://user-images.githubusercontent.com/10145585/50384454-c2b0ee80-06ea-11e9-949f-5d968278a4dd.png)|![boy_global_0128](https://user-images.githubusercontent.com/10145585/50384456-cb092980-06ea-11e9-9180-05bf4d08ecf8.png)|![boy_global_0188](https://user-images.githubusercontent.com/10145585/50384457-d2c8ce00-06ea-11e9-85b5-0e3f65d94b4b.png)|![boy_global_0254](https://user-images.githubusercontent.com/10145585/50384458-dd836300-06ea-11e9-83ed-24af0d8247b5.png)
|:---:|:---:|:---:|:---:|:---:|

**Quantitative Results**

OTB dataset

|![otb_precision_plot_scale](https://user-images.githubusercontent.com/10145585/50384498-6f8b6b80-06eb-11e9-8d4b-f474cce1f3ae.jpeg)|![otb_success_plot_scale](https://user-images.githubusercontent.com/10145585/50384508-8df16700-06eb-11e9-906e-0b4effc9da79.jpeg)|
|:---:|:---:|

VOT dataset

|![vot_precision_plot_scale](https://user-images.githubusercontent.com/10145585/50384510-a6fa1800-06eb-11e9-99e6-72712c7f367e.jpeg)|![vot_success_plot_scale](https://user-images.githubusercontent.com/10145585/50384512-aeb9bc80-06eb-11e9-9d20-f1ecb414306b.jpeg)|
|:---:|:---:|

# References

1. Zhen, Xiantong, et al. ”Multi-target regression via robust low-rank learning.” IEEE transactions on pattern analysis and machine intelligence 40.2 (2018): 497-504
2. Fhog -- https://github.com/pdollar/toolbox/blob/master/
3. OTB -- http://cvlab.hanyang.ac.kr/trackerbenchmark
4. VOT -- http://www.votchallenge.net/vot2017/
