/** 
 * @author Andrew Meyer  
 * @version 1.0  
 *  
 * Kalman Filter implementation for Rambling Rocket Club.  
 */

 template <class T, size_t N>
 struct Array {
    // Storage
    T data[N];

    static size_t length() { return N; }
    using type = T;

    // Item access
    T &operator[](size_t index) { return data[index]; }
    const T &operator[](size_t index) const { return data[index]; }

    // Iterators
    T *begin() { return &data[0]; }
    const T *begin() const { return &data[0]; }
    T *end() { return &data[N]; }
    const T *end() const { return &data[N]; }

    // Comparisons
    bool operator==(const Array<T, N> &rhs) const {
        if (this == &rhs)
            return true;
        for (size_t i = 0; i < N; i++)
            if ((*this)[i] != rhs[i])
                return false;
        return true;
    }
    bool operator!=(const Array<T, N> &rhs) const {
        return !(*this == rhs);
    }
};

class RocketFilter {
    static constexpr double t = 5.00; // modify this as needed
    static constexpr double a = 0.500; // Alpha value for kalman filter. Fine tune to values specific for rocket
    static constexpr double b = 0.400; // Beta value for kalman filter. Fine tune to values specific for rocket
    static constexpr double y = 0.100; // Gamma value for kalman filter. Fine tune to values specific for rocket
    double x_pos;
    double x_velo;
    double x_accel;
    double x_jerk;
    double y_pos;
    double y_velo;
    double y_accel;
    double y_jerk;
    double z_pos;
    double z_velo;
    double z_accel;
    double z_jerk;
    double measuredData;
    double measuredYData;
    double measuredZData;

    /**  
     * Constructor for RocketFilter.  
     * @param x_pos double representing the position in the x direction  
     * @param y_pos double representing the position in the y direction  
     * @param z_pos double representing the position in the z direction  
     * @param x_velo double representing the velocity in the x direction  
     * @param y_velo double representing the velocity in the y direction  
     * @param z_velo double representing the velocity in the z direction  
     * @param x_accel double representing the acceleration in x direction  
     * @param y_accel double representing the acceleration in y direction  
     * @param z_accel double representing the acceleration in z direction  
     * @param x_jerk double representing the jerk in x direction  
     * @param y_jerk double representing the jerk in y direction  
     * @param z_jerk double representing the jerk in z direction  
     */
public:RocketFilter(double x_pos, double y_pos, double z_pos, double x_velo, double y_velo, double z_velo,
                 double x_accel, double y_accel, double z_accel, double x_jerk, double y_jerk, double z_jerk,
                 double currXPos, double currYPos, double currZPos, double currXVelo, double currYVelo, 
                 double currZVelo, double currXAccel, double currYAccel, double currZAccel, double prevXPos,
                 double prevYPos, double prevZPos, double prevXVelo, double prevYVelo, double prevZVelo) {
                 this->x_pos = x_pos;
                 this->y_pos = y_pos;
                 this->z_pos = z_pos;
                 this->x_velo = x_velo;
                 this->y_velo = y_velo;
                 this->z_velo = z_velo;
                 this->x_accel = x_accel;
                 this->y_accel = y_accel;
                 this->z_accel = z_accel;
                 this->x_jerk = x_jerk;
                 this->y_jerk = y_jerk;
                 this->z_jerk = z_jerk;
                 this->currXPos = currXPos;
                 this->currYPos = currYPos;
                 this->currZPos = currZPos;
                 this->currXVelo = currXVelo;
                 this->currYVelo = currYVelo;
                 this->currZVelo = currZVelo;
                 this->currXAccel = currXAccel;
                 this->currYAccel = currYAccel;
                 this->currZAccel = currZAccel;
                 this->prevXPos = prevXPos;
                 this->prevYPos = prevYPos;
                 this->prevZPos = prevZPos;
                 this->prevXVelo = prevXVelo;
                 this->prevYVelo = prevYVelo;
                 this->prevZVelo = prevZVelo;
                 this->prevXAccel = 0;
                 this->prevYAccel = 0;
                 this->prevZAccel = 0;                 
                 this->measuredData = 0;
                 this->measuredYData = 0;
                 this->measuredZData = 0;
    }

    /*  
     * State extrapolation equations.  
     */

    /**  
     * Estimates x position of rocket after time t.  
     * @return double representing the estimated x position  
     */
private:double estimateXPos() {
        return estimatePos(this->x_pos, this->x_velo, this->x_accel);
    }
    /**  
     * Estimates y position of rocket after time t.  
     * @return double representing the estimated y position  
     */
private:double estimateYPos() {
        return estimatePos(this->y_pos, this->y_velo, this->y_accel);
    }
    /**  
     * Estimates z position of rocket after time t.  
     * @return double representing the estimated z position  
     */
private:double estimateZPos() {
        return estimatePos(this->z_pos, this->z_velo, this->z_accel);
    }

    /**  
    * Helper method for estimating position for each x, y, z.  
    * @param pos double representing position of the rocket in a dimension  
    * @param velo double representing velocity of the rocket in a dimension  
    * @param accel double representing acceleration of the rocket in a dimension  
    * @return updated state position of rocket after time t  
    */
private:double estimatePos(double pos, double velo, double accel) {
        return pos + velo * t + ((accel / 2) * t * t);
    }

    /**  
     * Estimates x velocity of the rocket after time t.  
     * @return double representing the x velocity after time t  
     */
private:double estimateXVelo() {
        return estimateVelo(this->x_velo, this->x_accel);
    }
    /**  
     * Estimates x velocity of the rocket after time t.  
     * @return double representing the x velocity after time t  
     */
private:double estimateYVelo() {
        return estimateVelo(this->y_velo, this->y_accel);
    }
    /**  
     * Estimates x velocity of the rocket after time t.  
     * @return double representing the x velocity after time t  
     */
private:double estimateZVelo() {
        return estimateVelo(this->z_velo, this->z_accel);
    }

    /**  
     * Helper method for estimating velocity for each x, y, z.  
     * @param velo double representing velocity of the rocket in a dimension  
     * @param accel double representing acceleration of the rocket in a dimension  
     * @return updated state velocity of rocket after time t  
     */
private:double estimateVelo(double velo, double accel) {
        return velo + accel * t;
    }

    // Note that the values for acceleration are treated as constant for now. This may change later

    /**  
     * Estimates x acceleration of the rocket after time t.  
     * @return double representing the the x acceleration after time t  
     */
private:double estimateXAccel() {
        return this->x_accel;
    }

    /** 
     * Estimates y acceleration of the rocket after time t. 
     * @return double representing the the y acceleration after time t 
     */
private:double estimateYAccel() {
        return this->y_accel;
    }

    /** 
     * Estimates z acceleration of the rocket after time t. 
     * @return double representing the the z acceleration after time t 
     */
private:double estimateZAccel() {
        return this->z_accel;
    }

    /*  
     * Retrieving data. Note that these methods will be implemented after everything else. Values will be measured  
     * arduino robot.  
     */

private:double getXPos() {
        return this->measuredData;
    }

private:double getYPos() {
        return this->measuredYData;
    }

private:double getZPos() {
        return this->measuredZData;
    }

    /*  
     * Update state equations.  
     */

    /**  
     * Returns an updated xPos for Rocket and updates xPos.  
     * @return Updated xPos for Rocket.  
     */
private:double updateXPos() {
        double data = (1 - a) * this->estimateXPos() + a * this->getXPos();
        return data;
    }
    /** 
     * Returns an updated yPos for Rocket and updates yPos. 
     * @return Updated yPos for Rocket. 
     */
private:double updateYPos() {
        double data = (1 - a) * this->estimateYPos() + a * this->getYPos();
        return data;
    }

    /** 
     * Returns an updated zPos for Rocket and updates zPos. 
     * @return Updated zPos for Rocket. 
     */
private:double updateZPos() {
        double data = (1 - a) * this->estimateZPos() + a * this->getZPos();
        return data;
    }

    /**  
     * Returns an updated xVelo for Rocket and updates xVelo.  
     * @return Updated xVelo for Rocket.  
     */
private:double updateXVelo() {
        double data = this->estimateXVelo() + b * ((this->getXPos() - this->estimateXPos()) / t);
        return data;
    }
    /** 
     * Returns an updated yVelo for Rocket and updates yVelo. 
     * @return Updated yVelo for Rocket. 
     */
private:double updateYVelo() {
        double data = this->estimateYVelo() + b * ((this->getYPos() - this->estimateYPos()) / t);
        return data;
    }
    /** 
     * Returns an updated zVelo for Rocket and updates zVelo. 
     * @return Updated zVelo for Rocket. 
     */
private:double updateZVelo() {
        double data = this->estimateZVelo() + b * ((this->getZPos() - this->estimateZPos()) / t);
        return data;
    }

    /**  
     * Returns an updated xAccel for Rocket and updates xAccel.  
     * @return Updated xAccel for Rocket  
     */
private:double updateXAccel() {
        double data = this->estimateXAccel() + y * ((this->getXPos() - this->estimateXPos()) / (0.5 * t * t));
        return data;
    }
    /** 
     * Returns an updated yAccel for Rocket and updates yAccel. 
     * @return updated yAccel for Rocket 
     */
private:double updateYAccel() {
        double data =  this->estimateYAccel() + y * ((this->getYPos() - this->estimateYPos()) / (0.5 * t * t));
        return data;
    }

    /** 
     * Returns an updated zAccel for Rocket and updates zAccel. 
     * @return Updated zAccel for Rocket 
     */
private:double updateZAccel() {
        double data = this->estimateZAccel() + y * ((this->getZPos() - this->estimateZPos()) / (0.5 * t * t));
        return data;
    }

    static double calculateNewVelo(double currAccel, double prevAccel, double prevVelo) {
        return prevVelo + ((currAccel + prevAccel) / 2) * t;
    }

    static double calculateNewPos(double currVelo, double prevVelo, double prevPos) {
        return prevPos + ((currVelo + prevVelo) / 2) * t;
    }

    Array<double, 9> arrayPrevs() {
        Array<double, 9> arr = {this->prevXPos, this->prevYPos, this->prevZPos,
                                this->prevXVelo, this->prevYVelo, this->prevZVelo,
                                this->prevXAccel, this->prevYAccel, this->prevZAccel};
        return arr;
    }

    Array<double, 9> arrayCurrs() {
        Array<double, 9> arr = {this->currXPos, this->currYPos, this->currZPos,
                                this->currXVelo, this->currYVelo, this->currZVelo,
                                this->currXAccel, this->currYAccel, this->currZAccel};
        return arr;
    }

public:
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    void updateAccels() {
      this->currXAccel = Serial.measure(); // find actual command later
      this->currYAccel = Serial.measure();
      this->currZAccel = Serial.measure();
    }

    void updateAllVals() {
        // if this doesnt work try taking an a double array of length 9 as a parameter
        //Indicies 6, 7, 8 are accel. 3, 4, 5 are velo. 0, 1, 2 are pos.
        //Array<double, 9> currs = this->arrayCurrs();
        //Array<double, 9> prevs = this->arrayPrevs();
        this->currXVelo = RocketFilter::calculateNewVelo(this->currXAccel, this->prevXAccel, this->prevXVelo);
        this->currYVelo = RocketFilter::calculateNewVelo(this->currYAccel, this->prevYAccel, this->prevYVelo);
        this->currZVelo = RocketFilter::calculateNewVelo(this->currZAccel, this->prevZAccel, this->prevZVelo);
        this->currXPos = RocketFilter::calculateNewPos(this->currXVelo, this->prevXVelo, this->prevXPos);
        this->currYPos = RocketFilter::calculateNewPos(this->currYVelo, this->prevYVelo, this->prevYPos);
        this->currZPos = RocketFilter::calculateNewPos(this->currZVelo, this->prevZVelo, this->prevZPos);
        this->prevXVelo = this->currXVelo;
        this->prevYVelo = this->currYVelo;
        this->prevZVelo = this->currZVelo;
        this->prevXAccel = this->currXAccel;
        this->prevYAccel = this->currYAccel;
        this->prevZAccel = this->currZAccel;
        // remove this last part if it breaks
        this->prevXPos = this->currXPos;
        this->prevYPos = this->currYPos;
        this->prevZPos = this->currZPos;

        // from here just use currPos to run filter
        this->measuredData = this->currXPos;
        this->measuredYData = this->currYPos;
        this->measuredZData = this->currZPos;
    }


    /** 
     * Setter for x_pos field. 
     * @param val value to be assigned to x_pos 
     */
    void setXPos(double val) {
        this->x_pos = val;
    }

    /** 
     * Setter for y_pos field. 
     * @param val value to be assigned to y_pos 
     */
    void setYPos(double val) {
        this->y_pos = val;
    }

    /** 
     * Setter for z_pos field. 
     * @param val value to be assigned to z_pos 
     */
    void setZPos(double val) {
        this->z_pos = val;
    }

    /** 
     * Setter for x_velo field. 
     * @param val value to be assigned to x_velo 
     */
    void setXVelo(double val) {
        this->x_velo = val;
    }

    /** 
     * Setter for y_velo field. 
     * @param val value to be assigned to y_velo 
     */
    void setYVelo(double val) {
        this->y_velo = val;
    }

    /** 
     * Setter for z_velo field. 
     * @param val value to be assigned to z_velo 
     */
    void setZVelo(double val) {
        this->z_velo = val;
    }

    /** 
     * Setter for x_accel field. 
     * @param val value to be assigned to x_accel 
     */
    void setXAccel(double val) {
        this->x_accel = val;
    }

    /** 
     * Setter for y_accel field. 
     * @param val value to be assigned to y_accel 
     */
    void setYAccel(double val) {
        this->y_accel = val;
    }

    /** 
     * Setter for z_accel field. 
     * @param val value to be assigned to z_accel 
     */
    void setZAccel(double val) {
        this->z_accel = val;
    }

    double getCurrXPos() {
        return this->currXPos;
    }

    double getCurrYPos() {
        return this->currYPos;
    }

    double getCurrZPos() {
        return this->currZPos;
    }

    /** 
     * Setter for measuredData field. 
     * @param val value to be assigned to measuredData 
     */
    void setMeasuredData(double val) {
        this->measuredData = val;
    }

    void setMeasuredYData(double val) {
        this->measuredYData = val;
    }

    void setMeasuredZData(double val) {
        this->measuredZData = val;
    }

    /**  
     * delivers an array of all of the estimated state fields of the rocket  
     * @return an array of length 9 containing all of the estimates state values for each dimension  
     */
    Array<double, 9> estimateState() {
      Array<double, 9> arr = {this->updateXPos(), this->updateYPos(), this->updateZPos(),
                         this->updateXVelo(), this->updateYVelo(), this->updateZVelo(),
                         this->updateXAccel(), this->updateYAccel(), this->updateZAccel()};
      return arr;
    }

    /**  
     * delivers an array of all of the predicted state fields of the rocket  
     * @return an array of length 9 containing all of the predicted state values for each dimension  
     */
    Array<double, 9> predictState() {
        Array<double, 9> arr = {this->updateXPos(), this->updateYPos(), this->updateZPos(),
                           this->updateXVelo(), this->updateYVelo(), this->updateZVelo(),
                           this->updateXAccel(), this->updateYAccel(), this->updateZAccel()};
        return arr;
    }

    /** 
     * Prints the data from a double array of length 9. 
     * @param arr array of length 9 which will have its contents printed 
     */
    static void printArray(Array<double, 9> arr) {
        for (double i: arr) {
            Serial.println(i);
        }
    }

    /** 
     * Test method for filter. 
     * @param i int representing the iteration to be tested 
     * @param measuredVals a double array containing the measured data values 
     */
    void filterTest(int i, Array<double, 10> measuredVals) {
        Serial.print("Printing iteration ");
        Serial.print(i);
        Serial.println(":");
        RocketFilter::printArray(this->estimateState());
        Serial.println("Printing updated state:");
        RocketFilter::printArray(this->predictState());
        Array<double, 9> states = this->predictState();

        //Updating values
        this->setXPos(states[0]);
        this->setYPos(states[1]);
        this->setZPos(states[2]);
        this->setXVelo(states[3]);
        this->setYVelo(states[4]);
        this->setZVelo(states[5]);
        this->setXAccel(states[6]);
        this->setYAccel(states[7]);
        this->setZAccel(states[8]);

        this->setMeasuredData(measuredVals[i]);

        Serial.print("iteration ");
        Serial.print(i);
        Serial.println(" finished.\n");
    }
};

RocketFilter* rocket = new RocketFilter(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
Serial.measure(), Serial.measure(), Serial.measure(), 0, 0, 0, 0, 0, 0);
double t = 5.00; // time between iterations
int i = 0;

  void setup() {
    Serial.println("Initial x, y, z pos: 0\nInitial x, y, z velo: 0\nInitial x, y, z accel: 0\n");
  }

  void loop() {
        Serial.print("Printing iteration ");
        Serial.print(i);
        Serial.println(":");

        rocket->updateAccels();
        rocket->updateAllVals();

        RocketFilter::printArray(rocket->estimateState());
        Serial.println("Printing updated state:");
        RocketFilter::printArray(rocket->predictState());
        Array<double, 9> states = rocket->predictState();

        //Updating values
        rocket->setXPos(states[0]);
        rocket->setYPos(states[1]);
        rocket->setZPos(states[2]);
        rocket->setXVelo(states[3]);
        rocket->setYVelo(states[4]);
        rocket->setZVelo(states[5]);
        rocket->setXAccel(states[6]);
        rocket->setYAccel(states[7]);
        rocket->setZAccel(states[8]);

        Serial.print("iteration ");
        Serial.print(i);
        Serial.println(" finished.\n");
        ++i;
        delay(t * 1000); // measurement in ms not s
  }

