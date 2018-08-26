# **CarND Term 2 Project 2 Writeup** 
# **Unscented Kalman Filter** 

#### This repository contains a c++ implementation of an unscented kalman filter. Please compare also to a previous repository with an [extended kalman filter](https://github.com/Anner-deJong/CarND-Extended-Kalman-Filter-Project/blob/master/writeup.md). The unscented filter is able to perform sensor fusion by combining incoming laser measurements and radar measurements,  The repository is structured so that it works together with the [simulation environment provided by Udacity](https://github.com/udacity/self-driving-car-sim/releases). (The link between the environment and the code is via [uWebSocketIO](https://github.com/uNetworking/uWebSockets), and already provided by Udacity). Furthermore, this writeup does *not* include any explanation about Kalman filters, only the current implementation details. A great Kalman filter theory resource can be found at [iLectureOnline](http://www.ilectureonline.com/lectures/subject/SPECIAL%20TOPICS/26/190).

---

## Important scripts

This writeup will give an overview of the extended kalman filter implementation by going through 3 important scripts:

#### [1. main.cpp](#1.-main.cpp),
#### [2. tools.cpp](#2.-tools.cpp),
#### [3. FusionEKF.cpp](#3.-FusionEKF.cpp), and
#### [4. kalman_filter.cpp](#4.-kalman_filter.cpp).

### 1. main.cpp

Exactly as with the extended kalman filter project, this script takes care of the uWebSocketIO server connection with the simulator, and passing measurements (in timesteps) from the simulator to the filter by calling `FusionEKF::ProcessMeasurement()` on an FusionEKF object. It also continuously updates the RMSE by taking all ground truths and predictions up to a certain timestep and passing them on to `Tools::CalculateRMSE()`.

### 2. tools.cpp

#### CalculateRMSE()

    VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,const vector<VectorXd> &ground_truth)

Upon finishing measurements and making a prediction at each timestep, the average Root Mean Square Error is calculated over all predictions and ground truths up to said timestep.
    
#### CalculateJacobian()

    MatrixXd Tools::CalculateJacobian(const VectorXd& x_state, const MatrixXd& prev_Hj_)

### 3. FusionEKF.cpp

This script implements a class (`FusionEKF`) that takes care of keeping around a Kalman filter object:
    
    KalmanFilter ekf_;
as well as updating it with new measurements:

    void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack)

Upon receiving the first measurement through `ProcessMeasurement()`, the measurement argument is used to initialize the Kalman filter object state (other filter matrices get initialized upon construction of the `FusionEKF` object itself, or get created during an update call).

In any following call to `ProcessMeasurement()` a measurement gets used to update the filter object. Again, depending on whether the measurement is from a radar or from a laser sensor, different measurement data is available and different update functions get called. A laser measurement is in cartesian coordinates, and a regular Kalman filter suffices: `KalmanFilter::Update()`. A radar measurement is in polar coordinates, and requires an Extended Kalman filter: `KalmanFilter::UpdateEKF()`.

### 4. kalman_filter.cpp

The kalman filter object has two purposes:

* Keeping the filter parameters around
* Updating the filter parameters with a newly received measurement

Updating is done based on the type of measurement received by the FusionEKF object, as stated above.
The update function for laser measurements (regular kalman filter) is:

    void KalmanFilter::Update(const Eigen::VectorXd &z);

The update function for radar measurements (extended kalman filter) is:

    void KalmanFilter::UpdateEKF(const Eigen::VectorXd &z);
    
Both of these functions only implement the parts of the update step which are different for the regular kalman filter as compared to the extended kalman filter. They both finish by internally calling the same function private `_UpdateCommon()` to fulfill the common part of the update step:
    
    void KalmanFilter::_UpdateCommon(const Eigen::VectorXd &y);
    
The `UpdateEKF()` function requires some extra care to prevent division by zero and non-normalized angles. For readability, angles are normalized with a separate private function, as suggested by a Udacity comment:

    void KalmanFilter::_NormalizeTanAngle(double& phi);

## Conclusion

Upon testing inside the simulator environment dataset 1, the RMSE output is well under the upper RMSE threshold set by Udacity:

* RMSE: [0.0973, 0.0855, 0.4513, 0.4399]
* Upper threshold: [.11, .11, 0.52, 0.52]

This means my first C++ project's implementation is correct and working well!

## Note

Currently, upon a restart or a dataset switch in the simulation environment, main.cpp does **not** forward any flag/event. This means that upon restart of the environment, the kalman filter object is **not** re-initialized. This results in the following behaviour:

* Not reinitializing them actually makes the performance worse on a second run as compared to a first run
* Not reinitializing upon switching environments can actually make the kalman filter code crash

I filed an issue in the original udacity repository asking if and how this would be possible.


