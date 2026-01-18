#!/usr/bin/env bash
set -euo pipefail

OUT=""
SAMPLE_SEC="1"

usage() {
  echo "Usage: monitor_cgroup.sh --out /path/to/stats.tsv [--sample_sec 1]" >&2
  exit 2
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --out) OUT="$2"; shift 2;;
    --sample_sec) SAMPLE_SEC="$2"; shift 2;;
    -h|--help) usage;;
    *) echo "Unknown arg: $1" >&2; usage;;
  esac
done

[[ -n "${OUT}" ]] || usage
SAMPLE_SEC="${SAMPLE_SEC//[[:space:]]/}"
if [[ -z "${SAMPLE_SEC}" ]]; then SAMPLE_SEC="1"; fi
if [[ "${SAMPLE_SEC}" == "0" ]]; then exit 0; fi

mkdir -p "$(dirname "${OUT}")"

CG="/sys/fs/cgroup"
CGV2=0
if [[ -f "${CG}/cgroup.controllers" ]]; then
  CGV2=1
fi

read_cpu_usage_usec() {
  if [[ "${CGV2}" -eq 1 && -r "${CG}/cpu.stat" ]]; then
    awk '$1=="usage_usec"{print $2; exit}' "${CG}/cpu.stat" 2>/dev/null || echo "0"
  elif [[ -r "${CG}/cpuacct/cpuacct.usage" ]]; then
    ns="$(cat "${CG}/cpuacct/cpuacct.usage" 2>/dev/null || echo 0)"
    awk -v ns="${ns}" 'BEGIN{printf "%.0f\n", ns/1000.0}'
  else
    echo "0"
  fi
}

read_mem_current() {
  if [[ "${CGV2}" -eq 1 && -r "${CG}/memory.current" ]]; then
    cat "${CG}/memory.current" 2>/dev/null || echo "0"
  elif [[ -r "${CG}/memory/memory.usage_in_bytes" ]]; then
    cat "${CG}/memory/memory.usage_in_bytes" 2>/dev/null || echo "0"
  else
    echo "0"
  fi
}

read_mem_peak() {
  if [[ "${CGV2}" -eq 1 && -r "${CG}/memory.peak" ]]; then
    cat "${CG}/memory.peak" 2>/dev/null || echo "0"
  elif [[ -r "${CG}/memory/memory.max_usage_in_bytes" ]]; then
    cat "${CG}/memory/memory.max_usage_in_bytes" 2>/dev/null || echo "0"
  else
    echo "0"
  fi
}

read_pids_current() {
  if [[ "${CGV2}" -eq 1 && -r "${CG}/pids.current" ]]; then
    cat "${CG}/pids.current" 2>/dev/null || echo "0"
  elif [[ -r "${CG}/pids/pids.current" ]]; then
    cat "${CG}/pids/pids.current" 2>/dev/null || echo "0"
  else
    echo "0"
  fi
}

echo -e "epoch_sec\tcpu_usage_usec\tmemory_current_bytes\tmemory_peak_bytes\tpids_current" > "${OUT}"

while true; do
  ts="$(date +%s)"
  cpu="$(read_cpu_usage_usec)"
  mem_cur="$(read_mem_current)"
  mem_peak="$(read_mem_peak)"
  pids="$(read_pids_current)"
  echo -e "${ts}\t${cpu}\t${mem_cur}\t${mem_peak}\t${pids}" >> "${OUT}"
  sleep "${SAMPLE_SEC}"
done
