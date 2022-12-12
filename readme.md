# IIR滤波器设计程序

树林错过了风儿，留下些树叶作为纪念。

## 依赖

1. GNU MP
2. GNU MPFR
3. MPREAL

## Build

```shell
mkdir build
cmake -B build . -G "Ninja"
cmake --build build
```

## Use
目前只写了椭圆滤波器,Butterworth、Chebyshev正在路上

滤波器类型:

| number | type     |
| ------ | -------- |
| 0      | Lowpass  |
| 1      | Highpass |
| 2      | Bandpass |
| 3      | Bandstop |

### 设计模拟滤波器

```
cd build
AFD_test.exe
```

### 设计数字滤波器

```
cd build
DFD_test.exe
```

输出中有系统的各项的LaTeX

带阻滤波器的设计和Mathematica设计出来的有较大出入(差异在于过渡带),Elliptic型似乎没有这个出入

