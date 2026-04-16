import time
import os
import psutil
import csv

# --------------------------------------------------
# Optional Streamlit support
# --------------------------------------------------
try:
    import streamlit as st
    STREAMLIT_AVAILABLE = True
except ImportError:
    STREAMLIT_AVAILABLE = False

# --------------------------------------------------
# External callback (API / GUI integration)
# set_step_callback(fn) donde fn(event, step_name, success)
# --------------------------------------------------
_step_callback = None

def set_step_callback(callback):
    """Registra un callback externo para recibir eventos de steps."""
    global _step_callback
    _step_callback = callback


def emit_qa_result(data: dict):
    """Emite el resultado estructurado de un QA al callback externo."""
    if _step_callback:
        try:
            _step_callback("qa_result", "", True, data)
        except Exception:
            pass

# --------------------------------------------------
# Paths (GUI-safe, CLI-safe)
# --------------------------------------------------
LOG_FILE = "user_data/logs/pipeline_ui.log"
CSV_FILE = "user_data/logs/pipeline_times.csv"


def _ensure_state():
    """
    Initialize state storage.
    Uses Streamlit session_state if available,
    otherwise falls back to a module-level dict.
    """
    if STREAMLIT_AVAILABLE:
        state = st.session_state
    else:
        if not hasattr(_ensure_state, "_state"):
            _ensure_state._state = {}
        state = _ensure_state._state

    state.setdefault("current_step", None)
    state.setdefault("current_step_start_time", None)
    state.setdefault("current_step_start", None)
    state.setdefault("pipeline_log", [])
    state.setdefault("pipeline_total_time", 0.0)

    return state


def start_step(step_name: str):
    """Mark the start of a pipeline step."""
    state = _ensure_state()

    state["current_step"] = step_name
    state["current_step_start"] = time.time()
    state["current_step_start_time"] = time.strftime(
        "%Y-%m-%d %H:%M:%S"
    )

    # Log buffer (GUI + CLI)
    state["pipeline_log"].append(f"▶ Started: {step_name}")
    state["pipeline_log"] = state["pipeline_log"][-50:]

    # External callback
    if _step_callback:
        try:
            _step_callback("start", step_name, True)
        except Exception:
            pass

    # Console output
    print(f"\n── {step_name}")


def end_step(success: bool = True):
    """Finish a step, log duration and system resources."""
    state = _ensure_state()

    step_name = state.get("current_step")
    start_epoch = state.get("current_step_start")
    start_str = state.get("current_step_start_time")

    if step_name is None or start_epoch is None:
        return

    # Duration
    elapsed = time.time() - start_epoch
    state["pipeline_total_time"] += elapsed

    # Status
    status_symbol = "✅" if success else "❌"
    msg = f"{status_symbol} {step_name} finished in {elapsed:.2f} s"
    state["pipeline_log"].append(msg)
    state["pipeline_log"] = state["pipeline_log"][-50:]

    # External callback
    if _step_callback:
        try:
            _step_callback("end", step_name, success)
        except Exception:
            pass

    # Console output
    if success:
        print(f"   ✅ {elapsed:.1f}s")
    else:
        print(f"   ❌ failed after {elapsed:.1f}s")

    # System resources
    ram = psutil.virtual_memory().used / (1024 ** 3)
    cpu = psutil.cpu_percent(interval=0.1)

    end_str = time.strftime("%Y-%m-%d %H:%M:%S")

    # --------------------------------------------------
    # 1. Log file (.log)
    # --------------------------------------------------
    try:
        os.makedirs(os.path.dirname(LOG_FILE), exist_ok=True)
        with open(LOG_FILE, "a", encoding="utf-8") as f:
            f.write(
                f"{start_str} | {end_str} | {status_symbol} {step_name} | "
                f"time={elapsed:.2f}s | RAM={ram:.2f}GB | CPU={cpu:.1f}%\n"
            )
    except Exception as e:
        if STREAMLIT_AVAILABLE:
            st.error(f"[LOG ERROR] {e}")
        else:
            print(f"[LOG ERROR] {e}")

    # --------------------------------------------------
    # 2. CSV file (.csv)
    # --------------------------------------------------
    try:
        os.makedirs(os.path.dirname(CSV_FILE), exist_ok=True)
        write_header = not os.path.exists(CSV_FILE)

        with open(CSV_FILE, "a", newline="", encoding="utf-8") as csvfile:
            writer = csv.writer(csvfile)

            if write_header:
                writer.writerow([
                    "start_time",
                    "end_time",
                    "step",
                    "duration_s",
                    "ram_gb",
                    "cpu_percent",
                ])

            writer.writerow([
                start_str,
                end_str,
                step_name,
                f"{elapsed:.2f}",
                f"{ram:.2f}",
                f"{cpu:.1f}",
            ])
    except Exception as e:
        if STREAMLIT_AVAILABLE:
            st.error(f"[CSV ERROR] {e}")
        else:
            print(f"[CSV ERROR] {e}")

    # Reset step
    state["current_step"] = None
    state["current_step_start"] = None
    state["current_step_start_time"] = None
