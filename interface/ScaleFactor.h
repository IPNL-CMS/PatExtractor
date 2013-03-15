#pragma once

class ScaleFactorCollection;

class ScaleFactor {

  public:
    ScaleFactor():
      m_value(1), m_error_low(0), m_error_high(0) { 
        m_array = new std::vector<double>(3);
      };

    ScaleFactor(double value, double error):
      m_value(value), m_error_low(error), m_error_high(error) {
        m_array = new std::vector<double>(3);
      };

    ScaleFactor(double value, double error_low, double error_high):
      m_value(value), m_error_low(error_low), m_error_high(error_high) {
        m_array = new std::vector<double>(3);
      };

    ScaleFactor(const ScaleFactor& other) {

      m_array = new std::vector<double>(3);

      m_value = other.m_value;
      m_error_low = other.m_error_low;
      m_error_high = other.m_error_high;
    }

    ScaleFactor& operator=(const ScaleFactor& other) {
      m_value = other.m_value;
      m_error_low = other.m_error_low;
      m_error_high = other.m_error_high;
      
      return *this;
    }

    ~ScaleFactor() {
      delete m_array;
    }

    double getValue() const {
      return m_value;
    }

    double getErrorLow() const {
      return m_error_low;
    }

    double getErrorHigh() const {
      return m_error_high;
    }

  private:

    friend class ScaleFactorCollection;

    ScaleFactor(const std::vector<double>& array) {
      m_value = array[0];
      m_error_low = array[1];
      m_error_high = array[2];
    }

    std::vector<double>* getBackingArray() const {
      (*m_array)[0] = m_value;
      (*m_array)[1] = m_error_low;
      (*m_array)[2] = m_error_high;

      return m_array;
    }

    double m_value;
    double m_error_low;
    double m_error_high;

    mutable std::vector<double>* m_array;
};

class ScaleFactorCollection {
  public:
    ScaleFactorCollection():
      m_owner(false), m_array(NULL) { };

    ~ScaleFactorCollection() {
      if (m_owner)
        delete m_array;

      m_array = NULL;
    }

    void push_back(const ScaleFactor& scaleFactor) {
      m_array->push_back(*scaleFactor.getBackingArray());
    }

    const ScaleFactor at(uint32_t index) const {
      return ScaleFactor((*m_array)[index]);
    }

    size_t size() const {
      return m_array->size();
    }

    std::vector<std::vector<double>>*& getBackingArray() {
      return m_array;
    }

    void clear() {
      m_array->clear();
    }

    void setWriteMode() {
      m_array = new std::vector<std::vector<double>>();
      m_owner = true;
    }

  private:
    bool m_owner;
    std::vector<std::vector<double>>* m_array;
};

