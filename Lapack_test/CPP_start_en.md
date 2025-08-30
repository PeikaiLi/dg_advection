# MSYS2 UCRT64 C++ Development Environment Guide (No mingw64)

Applicable scenarios:
- You want to develop modern C++ projects on Windows
- You only want to use **UCRT64**, and avoid mingw64 / msys
- You want to completely start over from scratch

---

##  Step 1: Uninstall Old MSYS2

1. Go to Control Panel → Uninstall MSYS2  
2. Delete the directory:  

   ```
   C:\msys64
   ```

---

## Step 2: Download and Install the Latest MSYS2

1. Open the official website: https://www.msys2.org  

2. Download the latest installer (example):  

   ```
   https://github.com/msys2/msys2-installer/releases/download/2025-06-22/msys2-x86_64-20250622.exe
   ```

3. Run the installer, keep the default path: `C:\msys64`

---

##  Step 3: Update MSYS2 Core System

1. Open:  

   ```
   MSYS2 MSYS
   ```

2. Enter the following command to update core packages:  

   ```
   pacman -Syu
   ```

3. If it asks you to close the window, just do so.  

4. Restart **MSYS2 MSYS**, then run:  

   ```
   pacman -Su
   ```

---

## Step 4: Enter the UCRT64 Environment

Open:  

```
MSYS2 MinGW UCRT64
```

You should see the terminal like this:  

```
Peikai_Li@YOUR_COMPUTER UCRT64 ~
```

Confirm and continue.

---

## Step 5: Install the UCRT64 C++ Toolchain

Run:  

```
pacman -S mingw-w64-ucrt-x86_64-toolchain
```

It will prompt you to select 13 components. Just press **Enter** to install all.

---

## Step 6: Verify Installation

1. Check g++ path:  

   ```
   which g++
   ```

   Expected:  

   ```
   /ucrt64/bin/g++
   ```

2. Check version:  

   ```
   g++ --version
   ```

3. Compile Hello World in `D:\document\practice_tem`  

   ```bash
   cd /D/document/practice_tem
   
   echo '#include <iostream>
   int main() { std::cout << "Hello UCRT64!" << std::endl; return 0; }' > hello.cpp
   
   g++ hello.cpp -o hello
   ./hello
   ```

   Expected output:  

   ```
   Hello
   ```

---

You can repeat install commands if needed:

```bash
pacman -Syu
pacman -S mingw-w64-ucrt-x86_64-toolchain
pacman -S mingw-w64-ucrt-x86_64-toolchain mingw-w64-ucrt-x86_64-openblas mingw-w64-ucrt-x86_64-lapack
pacman -S make
```

Terminal must be:  

```
MSYS2 UCRT64
```

g++ path must be:  

```
/ucrt64/bin/g++
```

Check installed packages:  

```bash
pacman -Q mingw-w64-ucrt-x86_64-openblas
pacman -Q mingw-w64-ucrt-x86_64-lapack
```

Example output:  

```
mingw-w64-ucrt-x86_64-openblas 0.3.30-2
mingw-w64-ucrt-x86_64-lapack 3.12.1-1
```

Make test:  

```bash
make clean
make
rm -f test.exe *.o
g++ -O2 -std=c++17 -o test.exe lapack_test.cpp -lopenblas
```

---

## Configure Windows System Environment Variables (UCRT64)

Add this path to Windows `Path` environment variable:  

```
C:\msys64\ucrt64\bin
```
