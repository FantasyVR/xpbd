# XPBD CPU

## 简介
这是一个简单的CPU版本的XPBD实现，目前仅使用了Distance Constraint和Volume Constraint

## 依赖
1. libigl(包含Eigen3， Embree3)

## 如何使用
1. `git clone http://oa.unidraw.com:9207/yupeng/xpbd_cpu.git`
2. `cd xpbd_cpu`
3. `cmake -S . -B build`
4. 打开 `build`文件夹下的 `xpbd_cpu.sln`开始编译

## 快捷键
`p` 仿真暂停/开始
`a`, `d`：x轴移动
`w`, `s`：y轴移动
`r`, `t`: z轴移动