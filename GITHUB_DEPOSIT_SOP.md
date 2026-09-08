# 手把手操作卡：还剩 3 件事（约 5 分钟）

更新时间：2026-09-08 19:2x

> 说明：GitHub 网页在你那边显示的是英文界面，所以下面每个按钮都写成
> **「英文原名（中文意思）」**，你照着找英文单词点就行。

## 先看：已经做完的（不用你动手）

| 项目 | 状态 |
|---|---|
| 分支 `main` 推送 | 71 个文件，最新提交 `4532024` |
| Tag `v1.0.0` | 已推送 |
| **GitHub Release v1.0.0** | 已自动创建，<https://github.com/yyx-4113/sleep-deprivation-scrna/releases/tag/v1.0.0> |
| 旧 `master` 分支 | 14.5 MB 旧内容原样保留，未删改 |
| 用户名纠正 | 全项目 `yongxinyang` → `yyx-4113` |

下面 3 件事**只能你本人做**（需要你的账号权限，我和本机网络都够不着），每步都写到"点哪个按钮"。

---

# 第 1 步：改仓库描述 ★最重要

> 现在仓库上写着 "…reveals Pomc as a central regulatory hub"（揭示了 Pomc 是中枢调控枢纽），
> 而你这篇稿子论证的正是"Pomc 中枢调控这个结论站不住"。
> 审稿人点开仓库会看到一句自我否定的话。**这一条优先级最高。**

### 1-1 打开设置页

浏览器地址栏粘贴这个网址、回车：

```
https://github.com/yyx-4113/sleep-deprivation-scrna/settings
```

**前提**：页面右上角要能看到你自己的头像。如果看到 **Sign in（登录）**，先登录再打开上面的链接。

### 1-2 找到描述输入框

进来后页面**最上方**会有一张卡片，标题是 **About（关于）**，里面有三项：

| 字段（英文原名） | 中文意思 | 现在的内容 |
|---|---|---|
| **Description** | 描述 | Single-cell transcriptomic analysis of sleep deprivation reveals Pomc as a central regulatory hub |
| Website | 网址 | （空） |
| Topics | 主题标签 | （空） |

> 小提示：Description 输入框**右侧有一个齿轮图标 ⚙**。
> 如果输入框是灰色不可编辑的，先点一下那个齿轮再改。

### 1-3 清空并粘贴新描述

把 Description 里**全部删掉**，然后粘贴这一行（一整行，不要换行）：

```
Reproducibility bundle for a data commentary on pseudoreplication, data leakage and inference limits in single-cell transcriptomics of sleep deprivation (GEO-driven inventory, code, adjudication records).
```

（这句的意思是：睡眠剥夺单细胞转录组研究中"伪重复、数据泄漏与推断局限"的数据评论文章之复现包。仓库描述面向国际读者，保留英文即可。）

### 1-4 保存

点输入框**右下角**的绿色按钮 **Save changes（保存更改）**。

### 1-5 怎么判断成功了

- 页面上方出现一条绿色提示条（写着 "Repository details updated"，即仓库信息已更新）
- 打开仓库主页 <https://github.com/yyx-4113/sleep-deprivation-scrna>，
  在仓库名 `sleep-deprivation-scrna` **正下方**，能看到刚粘贴的那句新描述

### 出问题了？

- **看不到 Settings（设置）标签** → 当前登录的账号不是 `yyx-4113`，或没登录。右上角头像 → **Sign out（退出）** → 重新用正确账号登录。
- **Save changes 点了没反应** → 检查描述是否超过 350 字符（上面这段约 190 字符，不会超）。
- **提示 You need admin access（需要管理员权限）** → 该仓库当时是用别的账号建的，换成有管理员权限的账号登录。

---

# 第 2 步：Zenodo 授权并拿 DOI

### 2-1 打开 Zenodo 并用 GitHub 登录

```
https://zenodo.org/
```

点页面**右上角**的 **Log in（登录）**，在登录方式里选带 GitHub 图标的按钮 **Log in with GitHub（用 GitHub 登录）**。

如果弹出授权页（标题 "Authorize Zenodo"，即授权 Zenodo），点绿色的 **Authorize zenodo（授权）**。

> ⚠️ 如果你这边 zenodo.org 打不开或一直转圈：直接跳到本页最后的 **「备选方案」**，
> 用 Science Data Bank（科学数据银行）一样能拿 DOI，不影响投稿。

### 2-2 打开 GitHub 仓库开关

登录后，点**右上角头像** → 菜单里点 **GitHub**。

或者直接走这个直达链接：

```
https://zenodo.org/account/settings/github/
```

页面会列出你 GitHub 账号下的所有仓库。找到 `yyx-4113/sleep-deprivation-scrna`：

- 列表长的话，用页面上的搜索框输入 `sleep-deprivation`
- 这一行**最右边有个开关**，点一下，让它从灰色变成 **ON（开启）**（蓝色/绿色）

### 2-3 等它归档（通常几分钟）

因为 **Release v1.0.0 已经建好了**，Zenodo 收到信号后会立刻开始归档。

想手动催一下也行：回到
<https://github.com/yyx-4113/sleep-deprivation-scrna/releases/tag/v1.0.0>
→ 右上角 **Edit（编辑）** → 什么都不改，直接拉到最下面点 **Update release（更新发布）**。

### 2-4 复制 DOI

几分钟后刷新 <https://zenodo.org/account/settings/github/>，找到这条记录点进去，
详情页右侧（或中部的 DOI 区块）会看到：

```
DOI 10.5281/zenodo.XXXXXXX
```

**复制这串号码发给我。**

### 2-5 一个容易踩的坑

详情页上如果有个开关叫 **"Cite all versions?"（引用所有版本？）**，**保持它是关闭状态**。
打开的话显示的是"跨版本 Concept DOI（概念 DOI）"，会随版本变动；
我们要的是锁定 `v1.0.0` 的 **Version DOI（版本 DOI）**，永久不变。

### 出问题了？

- **列表里找不到这个仓库** → 点页面上的 **sync now（立即同步）** 或刷新按钮；确认仓库是 Public（公开）。
- **一直没出现 DOI** → 确认 Release 已发布。打开
  <https://github.com/yyx-4113/sleep-deprivation-scrna/releases> 应能看到 v1.0.0，
  且标签不是 **Draft（草稿）**。

---

# 第 3 步：把默认分支从 master 改成 main（建议）

### 3-1 打开分支设置

```
https://github.com/yyx-4113/sleep-deprivation-scrna/settings/branches
```

### 3-2 切换

找到 **Default branch（默认分支）** 这一块。右侧有一个**双向箭头图标 ⇄**
（鼠标悬停会显示 "Switch default branch"，即切换默认分支）。

点它 → 在下拉框里选 **main** → 点 **Update（更新）**。

### 3-3 确认

弹出确认框，勾选/点 **I understand, update the default branch（我已知悉，更新默认分支）**。

### 3-4 怎么判断成功了

打开 <https://github.com/yyx-4113/sleep-deprivation-scrna>，
页面左上方分支下拉框显示的应该是 **main**（而不是 master），
文件列表里能看到 `CITATION.cff`、`LICENSE`、`README.md`、`scripts/`、`results/` 这些。

> 不用担心的：旧的 `master` 分支**没有被删除**，
> 在分支下拉框里仍可切换查看，2026-05 的旧内容完好无损。

---

# 改完怎么验收？（不用你自己判断）

你在网页改完后跟我说一声，我会跑这个脚本：

```bash
python scripts/check_repo_status.py
```

它通过 GitHub API（本机可达）逐项核对，输出**全中文**。例如刚刚跑出来的实况：

```
==============================================================
仓库状态核对：yyx-4113/sleep-deprivation-scrna
==============================================================
[未通过] 仓库描述
        当前：Single-cell transcriptomic analysis of sleep deprivation reveals Pomc as a central regulatory hub
[未通过] 默认分支 = master
[通过] 是否私有 = False
[通过] 现有分支 = ['main', 'master']
[通过] 最新 Release = v1.0.0（草稿=False）
[通过] main 分支含 CITATION.cff

==============================================================
仍有 2 项待办：
  - 第1步 仓库描述已更新
  - 第3步 默认分支为 main
```

全部变成「通过」就说明这一轮齐了。

---

# 拿到 DOI 之后，交给我

把 `10.5281/zenodo.XXXXXXX` 发我，我一次性做完：

1. 回填 `manuscript_commentary.md` 里的 `[DOI to be inserted]`（待填入 DOI）
2. 回填 `author_verification_statement.md`（作者核验声明）里的同一处占位
3. 重建 `manuscript_commentary.docx` 等全部 docx
4. 同步到仓库

---

# 备选方案：Zenodo 打不开时

用国内的 **Science Data Bank（科学数据银行）**，同样 mint 正式 DOI：

1. 打开 <https://www.scidb.cn> → 注册/登录
2. 新建数据集 → 填写标题、作者（Yang Yongxin）、摘要
3. 上传文件：先在浏览器打开下面这个链接下载打包好的源码（浏览器可访问 GitHub）

   ```
   https://github.com/yyx-4113/sleep-deprivation-scrna/archive/refs/tags/v1.0.0.tar.gz
   ```

   或者 zip 版：

   ```
   https://github.com/yyx-4113/sleep-deprivation-scrna/archive/refs/tags/v1.0.0.zip
   ```

4. 提交 → 审核通过后得到 DOI，发我即可

GigaDB（<https://gigadb.org>）也是同样流程，国际期刊认可度也高。

---

# 附：本机环境备忘（供以后复用）

- `https://github.com`(443) 与 `https://zenodo.org` 从本机**均不可达**；`api.github.com` 通；**SSH 22 可达** → git 一律走 SSH
- `gh` CLI 已装但未登录；本次推送与建 Release 都没依赖 gh 登录
- SSH 密钥：`~/.ssh/id_ed25519`
- Release 是用 GitHub Actions 自动建的（无 PAT 方案），见仓库 `.github/workflows/release.yml`
- 已验证做不到：Actions 的 token 改不了仓库 description / 默认分支（加 `administration: write` 也失败）
- `api.github.com` 偶发 504，核对脚本已内置 3 次重试
