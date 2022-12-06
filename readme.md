# 椭圆滤波器设计程序

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

## USE

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

