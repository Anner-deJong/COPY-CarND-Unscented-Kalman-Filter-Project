# **CarND Term 2 Project 2 Writeup** 
# **Unscented Kalman Filter** 

#### This repository contains a c++ implementation of an unscented kalman filter. Please compare also to a previous repository with an [extended kalman filter](https://github.com/Anner-deJong/CarND-Extended-Kalman-Filter-Project/blob/master/writeup.md). The unscented filter is able to perform sensor fusion by combining incoming laser measurements and radar measurements,  The repository is structured so that it works together with the [simulation environment provided by Udacity](https://github.com/udacity/self-driving-car-sim/releases). (The link between the environment and the code is via [uWebSocketIO](https://github.com/uNetworking/uWebSockets), and already provided by Udacity). Furthermore, this writeup does *not* include any explanation about Kalman filters, only the current implementation details. A great Kalman filter theory resource can be found at [iLectureOnline](http://www.ilectureonline.com/lectures/subject/SPECIAL%20TOPICS/26/190).

---

## Important scripts

This writeup will give an overview of the extended kalman filter implementation by going through 3 important scripts:

#### [1. main.cpp](#1.-main.cpp),
#### [2. tools.cpp](#2.-tools.cpp), and
#### [3. ukf.cpp](#3.-ukf.cpp).

### 1. main.cpp

Exactly as with the extended kalman filter project, this script takes care of the uWebSocketIO server connection with the simulator, and passing measurements (in timesteps) from the simulator to the filter by calling `ukf.ProcessMeasurement(meas_package)` on an FusionEKF object. It also continuously updates the RMSE by taking all ground truths and predictions up to a certain timestep and passing them on to `tools.CalculateRMSE(estimations, ground_truth);`.

### 2. tools.cpp

Upon finishing measurements and making a prediction at each timestep, the average Root Mean Square Error is calculated over all predictions and ground truths up to said timestep.

    VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,const vector<VectorXd> &ground_truth)

### 3. ukf.cpp

This script implements the unscented kalman filter class. The only public function of the class, `UKF::ProcessMeasurement()`, takes care of everything. On the first call, it will use the first measurement to initialize the state vector. Every subsequent call, it performs a full kalman filter cycle, i.e. a prediction step and an update step. The rest of the code is organized accordingly, and `UKF::ProcessMeasurement()` calls the correct functions for each of the two steps in chronological order.
    
    void UKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
    
        ////// Part A: PREDICTION STEP //////
        // Step 1: Generate (augmented) sigma points
        void UKF::AugmentSigmaPoints();

        // Step 2: Predict sigma points
        void UKF::PredictSigmaPoints();

        // Step 3: Predict state mean and covariance
        void UKF::PredictMeanAndCovariance();

        ////// Part B: UPDATE STEP //////
        // Step 1: Predict (radar/laser) measurements

        if (meas_pack.sensor_type_ == MeasurementPackage::RADAR) {
        void UKF::PredictRadarMeasurement(); // Radar measurement prediction
        }
        else if (meas_pack.sensor_type_ == MeasurementPackage::LASER) {
        void UKF::PredictLaserMeasurement(); // Lidar measurement prediction
        }

        // Step 2: Update state
        void UKF::UpdateState();
    }

As can be seen, except for *Part B Step 1*, all functions are the same for both radar as well as laser measurements. The difference in processing measurements from the two sensor modalities is handled by calling either `PredictRadarMeasurement()` or `PredictLaserMeasurement()`.

This script also holds all the information for covariance matrix standard deviations. Most of these values are provided by the  process noise parameters

## Conclusion

Upon testing inside the simulator environment dataset 1, the RMSE output is well underneath the RMSE values previously obtained through the extended kalman filter. Furthermore, the current RMSE values are also well beneath the threshold set by Udacity:

* (x loc , y loc , x speed, y_speed)
* [0.0669, 0.0821,  0.3269,  0.2313] - RMSE unscented kalman filter
* [0.0973, 0.0855,  0.4513,  0.4399] - RMSE extended kalman filter
* [0.0900, 0.1000,  0.4000,  0.3000] - Upper Udacity threshold


## Note

Results in the above paragraph *Conclusion* are all for Dataset 1 inside the simulator. For Dataset 2, the RMSE results are a bit worse [0.1057, 0.0630, 0.7182, 0.2882], especially the x speed, which jumps up at the beginning and has a hard time to recover. <br>
The fact the x direction speed RMSE jumps up at the beginning makes it seem like the model is somehow making a minus mistake (Dataset 2 starts out exactly as Dataset 1 yet in opposite x direction). The fact the filter can recover and that every loop runs the same exact code, makes this option less plausible however. <br>
Another option could be that the tweakable process noise parameters are too optimized for Dataset 1, and don't generalize that well to Dataset 2. The validity of this option is also weakened though as both Datasets seem to be mere mirrors of eachother, tracking an object with the same speed and curve characteristics (which would suggest same process noise parameters).

Although not a requirement for the project, I would like to try and calculate Normalized Innovation Squared (NIS) values in order to check the consistency of the chosen noise parameters. This is currently not added yet and would have to wait until I have more time.

Furthermore, same as with the extended kalman filter:
Currently, upon a restart or a dataset switch in the simulation environment, main.cpp does **not** forward any flag/event. This means that upon restart of the environment, the kalman filter object is **not** re-initialized. This results in the following behaviour:

* Not reinitializing them actually makes the performance worse on a second run as compared to a first run
* Not reinitializing upon switching environments can actually make the kalman filter code crash

I filed an issue in the original udacity repository asking if and how this would be possible.


