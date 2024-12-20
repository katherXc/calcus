#include <iostream>
#include <cmath>

class ADVar {
public:
    double val; // 值
    double der; // 導數

    // 預設建構子
    ADVar() : val(0.0), der(0.0) {}

    // 使用值和導數建構
    ADVar(double v, double d) : val(v), der(d) {}

    // 從 double 隱式轉換，只給值, 導數預設 0
    ADVar(double v) : val(v), der(0.0) {}

    // 運算子重載：加法
    ADVar operator+(const ADVar& other) const {
        return ADVar(val + other.val, der + other.der);
    }
    ADVar operator+(double c) const {
        return ADVar(val + c, der);
    }
    friend ADVar operator+(double c, const ADVar& x) {
        return ADVar(c + x.val, x.der);
    }

    // 運算子重載：減法
    ADVar operator-(const ADVar& other) const {
        return ADVar(val - other.val, der - other.der);
    }
    ADVar operator-(double c) const {
        return ADVar(val - c, der);
    }
    friend ADVar operator-(double c, const ADVar& x) {
        return ADVar(c - x.val, -x.der);
    }

    // 運算子重載：乘法
    ADVar operator*(const ADVar& other) const {
        // (f*g)' = f'*g + f*g'
        double newVal = val * other.val;
        double newDer = der * other.val + val * other.der;
        return ADVar(newVal, newDer);
    }
    ADVar operator*(double c) const {
        return ADVar(val * c, der * c);
    }
    friend ADVar operator*(double c, const ADVar& x) {
        return ADVar(c * x.val, c * x.der);
    }

    // 運算子重載：除法
    ADVar operator/(const ADVar& other) const {
        // (f/g)' = (f'*g - f*g') / g^2
        double newVal = val / other.val;
        double newDer = (der * other.val - val * other.der) / (other.val * other.val);
        return ADVar(newVal, newDer);
    }
    ADVar operator/(double c) const {
        return ADVar(val / c, der / c);
    }
    friend ADVar operator/(double c, const ADVar& x) {
        // (c/x)' = (0 - c*x') / x^2 = -c*x' / x^2
        return ADVar(c / x.val, (-c * x.der) / (x.val * x.val));
    }

    // 取負號
    ADVar operator-() const {
        return ADVar(-val, -der);
    }
};

// 數學函數擴充
ADVar sin(const ADVar& x) {
    return ADVar(std::sin(x.val), std::cos(x.val) * x.der);
}

ADVar cos(const ADVar& x) {
    return ADVar(std::cos(x.val), -std::sin(x.val) * x.der);
}

ADVar tan(const ADVar& x) {
    // tan'(x) = sec^2(x) = 1 + tan^2(x)
    double t = std::tan(x.val);
    return ADVar(t, (1.0 + t * t) * x.der);
}

ADVar asin(const ADVar& x) {
    // asin'(x) = 1 / sqrt(1 - x^2)
    return ADVar(std::asin(x.val), x.der / std::sqrt(1 - x.val * x.val));
}

ADVar acos(const ADVar& x) {
    // acos'(x) = -1 / sqrt(1 - x^2)
    return ADVar(std::acos(x.val), -x.der / std::sqrt(1 - x.val * x.val));
}

ADVar atan(const ADVar& x) {
    // atan'(x) = 1 / (1 + x^2)
    return ADVar(std::atan(x.val), x.der / (1 + x.val * x.val));
}

ADVar sinh(const ADVar& x) {
    // sinh'(x) = cosh(x)
    return ADVar(std::sinh(x.val), std::cosh(x.val) * x.der);
}

ADVar cosh(const ADVar& x) {
    // cosh'(x) = sinh(x)
    return ADVar(std::cosh(x.val), std::sinh(x.val) * x.der);
}

ADVar tanh(const ADVar& x) {
    // tanh'(x) = 1 - tanh^2(x)
    double th = std::tanh(x.val);
    return ADVar(th, (1.0 - th * th) * x.der);
}

ADVar exp(const ADVar& x) {
    return ADVar(std::exp(x.val), std::exp(x.val) * x.der);
}

ADVar log(const ADVar& x) {
    return ADVar(std::log(x.val), (1.0 / x.val) * x.der);
}

ADVar sqrt(const ADVar& x) {
    // sqrt'(x) = (1/(2*sqrt(x)))
    double s = std::sqrt(x.val);
    return ADVar(s, (1.0 / (2.0 * s)) * x.der);
}

// pow 函式： pow(f,g) = f^g
// 當 g 為常數時： (f^c)' = c * f^(c-1)*f'
ADVar pow(const ADVar& f, double c) {
    double newVal = std::pow(f.val, c);
    double newDer = c * std::pow(f.val, c - 1) * f.der;
    return ADVar(newVal, newDer);
}

// 當 g 也是 ADVar: (f^g)' = f^g * [ g'*ln(f) + g*(f'/f) ] (假設 f>0)
ADVar pow(const ADVar& f, const ADVar& g) {
    double newVal = std::pow(f.val, g.val);
    double newDer = newVal * (g.der * std::log(f.val) + g.val * (f.der / f.val));
    return ADVar(newVal, newDer);
}

int main() {
    std::cout << "這個程式會示範自動微分的功能。" << std::endl;
    std::cout << "請輸入一個數值 x，程式將計算 f(x) 及 f'(x)：\n";
    std::cout << "我們定義 f(x) = x^3 + sin(x)*exp(x) + log(x) + sqrt(x)\n";
    std::cout << "並自動計算在該 x 值下的函數值及其導數。" << std::endl;

    double input;
    std::cout << "請輸入 x 的值: ";
    std::cin >> input;

    // 定義自變數 x，導數為1 (表示對 x 微分)
    ADVar x(input, 1.0);

    // f(x) = x^3 + sin(x)*exp(x) + log(x) + sqrt(x)
    ADVar f = pow(x, 3) + sin(x) * exp(x) + log(x) + sqrt(x);

    std::cout << "f(" << input << ") = " << f.val << "\n";
    std::cout << "f'(" << input << ") = " << f.der << "\n";

    return 0;
}
