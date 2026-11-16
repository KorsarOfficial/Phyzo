# Phyzo

Локальный вычислительный инструмент для PDE на C++20. Windows x64, консоль и Win32 GUI.

Скалярный FEM, смешанная схема P2/P1 для несжимаемых течений в 2D/3D, косинусный метод Галёркина. Стационарные и нестационарные расчёты, проверка чувствительности к сетке, квадратуре и шагу времени. Архив расчёта сохраняет модель, параметры, сетку, диагностику, исходники и исполняемый файл для повторного запуска.

## Сборка

Нужны MSVC с поддержкой C++20, Windows SDK и CMake 3.24+. Зависимости устанавливаются отдельно; сборщик их не скачивает.

В каталоге `%PHYZO_SDK%\deps` должны находиться:

- `json.hpp` — nlohmann/json 3.12.0;
- `eigen-3.4.0\Eigen` и `boost_1_86_0\boost` — Eigen 3.4.0 и Boost 1.86.0;
- `sundials\include` и `sundials\lib` — SUNDIALS 7.4.0, статическая Release-сборка x64 с CVODE, serial NVector, dense SUNMatrix и dense SUNLinearSolver;
- `miniz` — исходники miniz 3.0.2; `tinyxml2` — исходники tinyxml2 11.0.0.

Из корня репозитория в x64 Native Tools Command Prompt:

```bat
set PHYZO_SDK=C:\Phyzo-native
cl /std:c++20 /EHsc /utf-8 build.cpp /Fe:BuildPhyzo.exe
BuildPhyzo.exe
```

Результат: `%PHYZO_SDK%\bin\Phyzo.exe`. Промежуточные файлы сборки создаются в SDK.

## Запуск

```bat
Phyzo.exe examples
Phyzo.exe run heat3
Phyzo.exe stat nseq2
Phyzo.exe study heat3 --kind mesh --levels 3
Phyzo.exe selftest --count 1000 --seed 20261115 -o tests.json
Phyzo.exe --help
```

Вместо имени встроенного примера можно передать путь к JSON-модели. Без аргументов открывается GUI. `PHYZO_DATA` задаёт каталог результатов; `verify` проверяет архив, `replay` повторяет расчёт.
