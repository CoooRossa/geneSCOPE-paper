#!/usr/bin/env bash
set -euo pipefail

OUT=""
SAMPLE_SEC="1"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --out) OUT="$2"; shift 2;;
    --sample_sec) SAMPLE_SEC="$2"; shift 2;;
    *) echo "Unknown arg: $1" >&2; exit 1;;
  esac
done

if [[ -z "${OUT}" ]]; then
  echo "Usage: monitor_cgroup.sh --out <path.tsv> --sample_sec <int>" >&2
  exit 2
fi

SAMPLE_SEC="$(printf "%s" "${SAMPLE_SEC}" | tr -d " ")"
if [[ -z "${SAMPLE_SEC}" ]]; then SAMPLE_SEC="1"; fi
if [[ "${SAMPLE_SEC}" == "0" ]]; then exit 0; fi

mkdir -p "$(dirname "${OUT}")"

CGV2=0
if [[ -f /sys/fs/cgroup/cgroup.controllers ]]; then
  CGV2=1
fi

read_mem_current() {
  if [[ "${CGV2}" -eq 1 && -r /sys/fs/cgroup/memory.current ]]; then
    cat /sys/fs/cgroup/memory.current
  elif [[ -r /sys/fs/cgroup/memory/memory.usage_in_bytes ]]; then
    cat /sys/fs/cgroup/memory/memory.usage_in_bytes
  else
    echo ""
  fi
}

read_mem_peak() {
  if [[ "${CGV2}" -eq 1 && -r /sys/fs/cgroup/memory.peak ]]; then
    cat /sys/fs/cgroup/memory.peak
  elif [[ -r /sys/fs/cgroup/memory/memory.max_usage_in_bytes ]]; then
    cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes
  else
    echo ""
  fi
}

read_cpu_usage_usec() {
  if [[ "${CGV2}" -eq 1 && -r /sys/fs/cgroup/cpu.stat ]]; then
    awk '$1=="usage_usec"{print $2}' /sys/fs/cgroup/cpu.stat
  elif [[ -r /sys/fs/cgroup/cpuacct/cpuacct.usage ]]; then
    ns="$(cat /sys/fs/cgroup/cpuacct/cpuacct.usage)"
    awk -v ns="${ns}" 'BEGIN{printf "%.0f\n", ns/1000.0}'
  else
    echo ""
  fi
}

if [[ ! -f "${OUT}" ]]; then
  echo -e "epoch\tiso_time\tmem_current_bytes\tmem_peak_bytes\tcpu_usage_usec\tcpu_usage_delta_usec" > "${OUT}"
fi

prev_cpu="$(read_cpu_usage_usec || true)"
if [[ -z "${prev_cpu}" ]]; then prev_cpu=""; fi

while true; do
  epoch="$(date +%s)"
  iso="$(date -Is)"
  mem_cur="$(read_mem_current || true)"
  mem_peak="$(read_mem_peak || true)"
  cpu_now="$(read_cpu_usage_usec || true)"

  delta=""
  if [[ -n "${cpu_now}" && -n "${prev_cpu}" && "${cpu_now}" =~ ^[0-9]+$ && "${prev_cpu}" =~ ^[0-9]+$ ]]; then
    if (( cpu_now >= prev_cpu )); then
      delta="$((cpu_now - prev_cpu))"
    fi
  fi
  prev_cpu="${cpu_now}"

  echo -e "${epoch}\t${iso}\t${mem_cur}\t${mem_peak}\t${cpu_now}\t${delta}" >> "${OUT}"
  sleep "${SAMPLE_SEC}"
done
