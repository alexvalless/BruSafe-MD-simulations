# Paso a paso para cada computadora del laboratorio

Todos los comandos van en la terminal **MSYS2 UCRT64** (menú Inicio → "MSYS2 UCRT64"),
**de uno en uno**, tal cual, sin cambiar nada (salvo el número de máquina del paso 8).
Primera vez en una PC: unos 25 minutos. Después, el trabajo corre solo (1 a 2 horas).

Lo que ya está resuelto y no hay que tocar: la versión de GROMACS es **2026.3 con CUDA**
(la misma en todas), la ruta de OneDrive está guardada en los scripts y los resultados se
copian solos al terminar.

---

## 1. GPU y Git

```bash
git --version; nvidia-smi --query-gpu=name --format=csv,noheader
```

Debe salir la versión de Git y `NVIDIA GeForce RTX 4070`. Si la GPU no sale, esa PC no sirve.

## 2. Instalar Python (sin administrador)

```bash
winget.exe install -e --id Python.Python.3.12 --scope user --accept-package-agreements --accept-source-agreements
```

Si `winget` no existe o falla: instalador de python.org (3.12, 64 bits), **"Install for me only"**
y "Add python.exe to PATH" marcado.

**Cierra la terminal y ábrela de nuevo.** Luego:

```bash
python --version
```

Debe decir `Python 3.12.x`. Usa `python`, no `python3` (ese es un atajo de la Microsoft Store
que no sirve; los scripts del proyecto ya lo ignoran).

## 3. Instalar numpy

```bash
python -m pip install --user numpy
```
```bash
python -c "import numpy; print(numpy.__version__)"
```

Debe imprimir una versión (por ejemplo `2.5.3`).

## 4. Clonar el repositorio

```bash
cd ~ && git clone https://github.com/alexvalless/BruSafe-MD-simulations.git && cd BruSafe-MD-simulations && git checkout claude/gallant-babbage-ydyxr2
```

## 5. El zip de GROMACS

```bash
ls -lh "/c/Users/A01563079/OneDrive - Instituto Tecnologico y de Estudios Superiores de Monterrey/BruSafe-MD/software/"
```

Debe aparecer `gromacs-2026.3-win64-cuda12.6.3-sm89.zip` de **unos 242 MB**.

- **Aparece pero no muestra el tamaño, o pesa 0:** OneDrive solo tiene el nombre, no el archivo.
  En el Explorador de Windows: clic derecho sobre el zip → **"Mantener siempre en este
  dispositivo"** y espera la marca verde sólida.
- **No existe la carpeta:** esa PC no tiene OneDrive iniciado. Baja el zip a `Descargas` con tu
  sesión de GitHub iniciada (y ciérrala al terminar):
  https://github.com/alexvalless/BruSafe-MD-simulations/actions/runs/37343539360/artifacts/11363541517
- Comprobación de integridad (opcional): la suma debe ser
  `ed5fd8e86a8733433a47cba8f816955070d20369ef10fbff6e77a6f38760a940`.

## 6. Preparar y verificar la máquina

```bash
./scripts/setup_machine.sh
```

Extrae GROMACS, deja bien tu `~/.bashrc` y revisa todo. Debe terminar con
`READY: CHI-12104-PC-XX`. Si dice `NOT READY`, te lista lo que falta.
**Cierra la terminal y abre una nueva** antes de seguir.

## 7. Revisar el plan

Los trabajos están en `docs/PLAN.txt`. Hoy son 7:

| Máquina | Trabajo | Duración aprox. |
|---|---|---|
| 0 | S8 (RNA libre) réplica 1 | 1 h |
| 1 | S8 réplica 2 | 1 h |
| 2 | S8 réplica 3 | 1 h |
| 3 | S4 (complejo) réplica 2 | 1.4 h |
| 4 | S4 réplica 3 | 1.4 h |
| 5 | S4 réplica 4 | 1.4 h |
| 6 | S4 réplica 5 | 1.4 h |

Apunta qué PC toma qué número para no repetir trabajos. (S4 réplica 1 ya está hecha.)

## 8. Lanzar (Terminal 1)

Copia **solo la línea de tu número**. Cada una actualiza el repo y lanza:

```bash
cd ~/BruSafe-MD-simulations && git pull --no-rebase && ./scripts/run_slot.sh 0 16 > slot_0.log 2>&1
```
```bash
cd ~/BruSafe-MD-simulations && git pull --no-rebase && ./scripts/run_slot.sh 1 16 > slot_1.log 2>&1
```
```bash
cd ~/BruSafe-MD-simulations && git pull --no-rebase && ./scripts/run_slot.sh 2 16 > slot_2.log 2>&1
```
```bash
cd ~/BruSafe-MD-simulations && git pull --no-rebase && ./scripts/run_slot.sh 3 16 > slot_3.log 2>&1
```
```bash
cd ~/BruSafe-MD-simulations && git pull --no-rebase && ./scripts/run_slot.sh 4 16 > slot_4.log 2>&1
```
```bash
cd ~/BruSafe-MD-simulations && git pull --no-rebase && ./scripts/run_slot.sh 5 16 > slot_5.log 2>&1
```
```bash
cd ~/BruSafe-MD-simulations && git pull --no-rebase && ./scripts/run_slot.sh 6 16 > slot_6.log 2>&1
```

Esa terminal queda **ocupada y en silencio** hasta terminar.

## 9. Ver el avance (Terminal 2, una ventana UCRT64 aparte)

El mismo comando en todas las máquinas; se refresca solo cada 30 segundos:

```bash
cd ~/BruSafe-MD-simulations && ./scripts/progress.sh watch
```

| Etapa | Cuánto | Qué ves |
|---|---|---|
| construcción + equilibración | ≈10 min | etapa `equil`; en `slot_N.log`: `building`, `equilibration` |
| producción | S8 ≈1 h, S4 ≈1.3 h | `prod`, ns del último checkpoint (cada 15 min), tamaño de `prod.xtc` y hace cuánto creció |
| postproceso + empaquetado | unos minutos | `post-processing`, `results packed` |

Salir de la Terminal 2 con **Ctrl+C** es seguro: no afecta el trabajo.

## 10. Reglas mientras corre

1. **Nunca Ctrl+C en la Terminal 1.** (Es lo que detuvo la PC-02 la primera vez.)
2. No cierres esa ventana.
3. Bloquea con **Win+L**; no cierres sesión y evita que el equipo se suspenda.
4. Si se corta la luz o se reinicia: vuelve a abrir UCRT64 y repite **el mismo comando del
   paso 8**. Retoma desde el último checkpoint y no repite lo ya terminado.

## 11. Al terminar

El prompt de la Terminal 1 regresa. Comprueba:

```bash
tail -3 ~/BruSafe-MD-simulations/slot_N.log
```

(con tu número en lugar de `N`). Debe decir `all jobs of slot N completed`.
Los resultados quedan en tu OneDrive, en `BruSafe-MD/<nombre-de-la-PC>`. Espera a que sus
archivos tengan la **marca verde** antes de cerrar sesión, y cierra la sesión de OneDrive y
de GitHub en el navegador (son computadoras compartidas).

---

## Si algo falla

| Síntoma | Qué hacer |
|---|---|
| `missing required command: gmx` | `git pull --no-rebase` (ya lo corrige) y relanza; o `source ~/.bashrc` |
| `Received the INT signal` en el log | Alguien pulsó Ctrl+C. Relanza el mismo comando del paso 8 |
| `NOT READY` | Lee la línea `[FAIL]` y corrígela; repite `setup_machine.sh` |
| `python3: Python was not found` | Usa `python` (paso 2) |
| El zip no muestra tamaño | Paso 5: "Mantener siempre en este dispositivo" |
| La trayectoria no crece en >10 min | Avisa con el número de máquina |
| `NO CHARMM-GUI build for ...` | Ese sistema aún no está listo; no es un error de tu PC |

Para pedir ayuda, copia: el nombre de la PC, tu número de máquina y la salida de

```bash
grep -E "^\[|Fatal|FATAL|rror" ~/BruSafe-MD-simulations/slot_N.log | tail -15
```
