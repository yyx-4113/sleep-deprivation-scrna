#!/usr/bin/env python3
"""核对复现仓库的公开状态（输出为中文）。

本机 github.com(443) 不可达，但 api.github.com 可达，故用 API 核对。
作者在网页端改完后运行本脚本即可逐项验收。

用法：  python scripts/check_repo_status.py
"""

import json
import os
import ssl
import sys
import time
import urllib.request

TOKEN = os.environ.get("GH_TOKEN")

REPO = "yyx-4113/sleep-deprivation-scrna"
API = f"https://api.github.com/repos/{REPO}"
WANT_DESC = (
    "Reproducibility bundle for a data commentary on pseudoreplication, data leakage "
    "and inference limits in single-cell transcriptomics of sleep deprivation "
    "(GEO-driven inventory, code, adjudication records)."
)

CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE


def get(url, retries=3):
    """带重试的请求：api.github.com 偶发 504，重试通常即可恢复。"""
    last = None
    for i in range(retries):
        try:
            headers = {"User-Agent": "repo-status-check"}
            if TOKEN:
                headers["Authorization"] = f"Bearer {TOKEN}"
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, context=CTX, timeout=30) as r:
                return json.load(r)
        except Exception as e:  # noqa: BLE001
            last = e
            if i < retries - 1:
                time.sleep(3)
    raise last


def mark(flag):
    return "通过" if flag else "未通过"


def main():
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:  # noqa: BLE001
        pass

    try:
        repo = get(API)
    except Exception as e:  # noqa: BLE001
        print(f"无法连接 api.github.com：{e}")
        return 1

    results = []

    print("=" * 62)
    print(f"仓库状态核对：{REPO}")
    print("=" * 62)

    # 1. 描述
    desc = repo.get("description") or ""
    same = desc.strip() == WANT_DESC.strip()
    results.append(("第1步 仓库描述已更新", same))
    print(f"[{mark(same)}] 仓库描述")
    print(f"        当前：{desc or '(空)'}")
    if not same:
        print(f"        应为：{WANT_DESC}")

    # 2. 默认分支
    dflt = repo.get("default_branch")
    is_main = dflt == "main"
    results.append(("第3步 默认分支为 main", is_main))
    print(f"[{mark(is_main)}] 默认分支 = {dflt}")

    # 3. 公开性
    priv = repo.get("private")
    results.append(("仓库为公开", priv is False))
    print(f"[{mark(priv is False)}] 是否私有 = {priv}")

    # 4. 分支完整性
    branches = [b["name"] for b in get(API + "/branches")]
    has_both = "main" in branches and "master" in branches
    results.append(("main 与 master 并存（旧内容未删）", has_both))
    print(f"[{mark(has_both)}] 现有分支 = {branches}")

    # 5. Release（用 release 列表而非 /latest，因为 /latest 会跳过 prerelease）
    try:
        rels = get(API + "/releases")
        rel = next((r for r in rels if r.get("tag_name") == "v1.0.0"), None)
        if rel is None:
            raise KeyError("v1.0.0 不在 release 列表中")
        good = not rel.get("draft")
        results.append(("Release v1.0.0 已发布", good))
        print(f"[{mark(good)}] Release v1.0.0（草稿={rel.get('draft')}）")
        print(f"        {rel.get('html_url')}")
    except Exception:  # noqa: BLE001
        results.append(("Release v1.0.0 已发布", False))
        print("[未通过] 未找到 Release")

    # 6. CITATION.cff
    try:
        get(f"{API}/contents/CITATION.cff?ref=main")
        results.append(("main 分支含 CITATION.cff", True))
        print("[通过] main 分支含 CITATION.cff")
    except Exception:  # noqa: BLE001
        results.append(("main 分支含 CITATION.cff", False))
        print("[未通过] main 分支缺少 CITATION.cff")

    failed = [n for n, v in results if not v]
    print("\n" + "=" * 62)
    if failed:
        print(f"仍有 {len(failed)} 项待办：")
        for n in failed:
            print("  -", n)
    else:
        print("全部通过，本轮操作已齐备。")
    return 0 if not failed else 2


if __name__ == "__main__":
    sys.exit(main())
