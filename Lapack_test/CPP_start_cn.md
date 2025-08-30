# MSYS2 UCRT64 C++ 开发环境完整重装指南（纯净、无 mingw64）

适用场景：
- 你想用 Windows 写现代 C++ 项目
- 你只用 UCRT64，不想碰 mingw64 / msys
- 你想完全从零重新来一次

---

## 🧹 第一步：卸载旧版 MSYS2

1. 控制面板 → 卸载 MSYS2
2. 删除目录：
   
        C:\msys64

---

## 📥 第二步：下载并安装最新 MSYS2

1. 打开官网：https://www.msys2.org

2. 下载最新版安装器（例如）：

        https://github.com/msys2/msys2-installer/releases/download/2025-06-22/msys2-x86_64-20250622.exe

3. 双击安装，路径保持默认：`C:\msys64`

---

## 🔄 第三步：更新 MSYS2 核心系统

1. 打开：

        MSYS2 MSYS

2. 输入以下命令更新核心包：

        pacman -Syu

3. 如果它提示你要关闭窗口，就照做。

4. 重启 MSYS2 MSYS，继续输入：

        pacman -Su

---

## ✅ 第四步：进入 UCRT64 环境

打开：

    MSYS2 MinGW UCRT64

你会看到终端变成这样：

    Peikai_Li@你的电脑名 UCRT64 ~

确认无误后进入下一步。

---

## 📦 第五步：安装 UCRT64 的 C++ 工具链

执行以下命令：

    pacman -S mingw-w64-ucrt-x86_64-toolchain

系统会提示你选择 13 个组件，直接 **按 Enter 全部安装**，不用输入数字。

---

## 🧪 第六步：测试是否安装成功

1. 查看 g++ 安装路径：

        which g++

    应该输出：

        /ucrt64/bin/g++

2. 查看版本：

        g++ --version

3. 编译一个 Hello World in D:\document\practice_tem

    `cd /D/document/practice_tem`

        echo '#include <iostream>
        int main() { std::cout << "Hello UCRT64!" << std::endl; return 0; }' > hello.cpp
        
        g++ hello.cpp -o hello
        
        ./hello

    你应该看到：

        Hello

---

你可以重复安装命令：

```bash
pacman -Syu
pacman -S mingw-w64-ucrt-x86_64-toolchain
pacman -S mingw-w64-ucrt-x86_64-toolchain mingw-w64-ucrt-x86_64-openblas mingw-w64-ucrt-x86_64-lapack
pacman -S make
```




终端必须是：

    MSYS2 UCRT64

路径必须是：

    /ucrt64/bin/g++

```bash
pacman -Q mingw-w64-ucrt-x86_64-openblas
pacman -Q mingw-w64-ucrt-x86_64-lapack

# check 
Peikai_Li@LPK UCRT64 ~
$ pacman -Q mingw-w64-ucrt-x86_64-openblas

mingw-w64-ucrt-x86_64-openblas 0.3.30-2

Peikai_Li@LPK UCRT64 ~
$

Peikai_Li@LPK UCRT64 ~
$ pacman -Q mingw-w64-ucrt-x86_64-lapack

mingw-w64-ucrt-x86_64-lapack 3.12.1-1


Peikai_Li@LPK UCRT64 /d/document/cpp_leaning/cpp_tutorior/Lapack_test
$ make clean
make
rm -f test.exe *.o
g++ -O2 -std=c++17 -o test.exe lapack_test.cpp -lopenblas

```



## ✅ 配置 Windows 系统环境变量（UCRT64）

将以下路径添加到系统环境变量 Path 中：`C:\msys64\ucrt64\bin`

