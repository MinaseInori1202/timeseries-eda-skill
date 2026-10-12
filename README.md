# timeseries-eda

時間序列探索性分析（EDA）與建模銜接的 **Claude Skill**，以 R 語言撰寫，涵蓋傳統統計與現代 tidyverts 兩條路線。

給 Claude 一句話說明資料與目的，它會先**確認你的分析取向**，再依對應路線走完整套流程，
最後交付一份可重跑的 Quarto 報告（`.qmd` 原始檔 + 編譯好的 `.html`）。

---

## 特色

- **開場先問取向，不自作主張** —— 模型的可解釋性 vs 預測準確／pipeline 穩健，兩條路線的工具與判讀完全不同
- **兩條完整路線** —— 傳統 Box–Jenkins（ADF、ACF/PACF/EACF、`forecast`）與現代 tidyverts（`tsibble`/`feasts`/`fable`）
- **金融序列專用分支** —— 股價一律轉對數報酬、六項摘要統計、ARCH effect 檢定、GARCH 方向判斷
- **殘差兩層診斷** —— 一階矩（Ljung–Box）與二階矩（ARCH），**兩層結論都必須寫進報告**
- **長季節不踩雷** —— 日資料年週期不套 SARIMA(365)，改用 STL 或動態諧波迴歸（Fourier）
- **季節差分有檢定依據** —— 月、季資料另做季節單根檢定決定 D，不只靠 ADF 或 `nsdiffs()` 預設值
- **外生變數（ARIMAX）** —— 迴歸誤差項為 ARMA，附三步驟檢查（迴歸 → 殘差單根／序列相關 → 聯合估計）與 R 套件對照
- **檢定寫法有規範** —— Ljung–Box 必須寫出 H₀／H₁、滯後期數與結果（依 Tsay 講義：一般 m = 10；GARCH 後標準化殘差及其平方各報 Q(10)、Q(20)）
- **報告格式已規範** —— 三線表、中英分開字型、單一自足 HTML、預測圖必附圖例、圖內不放主標題、預設不顯示程式碼

---

## 環境需求

| 項目 | 說明 |
|---|---|
| R | 4.x 以上 |
| Quarto | 用來把 `.qmd` 編譯成 `.html`（[下載](https://quarto.org/docs/get-started/)）|

R 套件（依實際路線安裝，不必全裝）：

```r
# 共通
install.packages(c("flextable", "e1071", "gridExtra"))

# 傳統路線
install.packages(c("forecast", "tseries", "TSA", "zoo", "urca", "lmtest", "sandwich"))

# 現代路線
install.packages("fpp3")          # 內含 tsibble / feasts / fable

# 金融序列（股價、報酬、波動度）
install.packages(c("quantmod", "FinTS", "rugarch", "fGarch"))

# 選配
install.packages(c("prophet", "gtsummary", "imputeTS", "uroot"))   # uroot：HEGY 季節單根檢定（尚未實測）
```

> **注意**：`library(fpp3)` 之後再 `library(FinTS)` 會遮蔽 `fable` 的 `ARIMA()`。
> 請**先載入 `FinTS`**，或改用命名空間呼叫 `FinTS::ArchTest()`。

---

## 安裝

### 方法一：請 Claude 幫你裝（最簡單）

1. 點本頁綠色 **Code** → **Download ZIP**
2. 開啟 Claude，直接跟它說：

```
我下載了 timeseries-eda-skill-main.zip，請幫我安裝這個 skill
```

Claude 會自動解壓並放到正確位置。

### 方法二：手動放置

解壓後把資料夾放到 Claude 的 skills 目錄，並命名為 `timeseries-eda`：

```
~/.claude/skills/timeseries-eda/     # 個人層級（所有專案都可用）
.claude/skills/timeseries-eda/       # 專案層級（只在該專案可用）
```

Windows 的 `~` 是 `C:\Users\<你的帳號>`。

### 確認安裝成功

資料夾結構應如下，且 **`SKILL.md` 必須在最外層**：

```
timeseries-eda/
├── SKILL.md
├── assets/
└── references/
```

放好後**重新啟動 Claude**，skill 才會生效。

---

## 使用方式

直接用自然語言說明**資料**與**目的**：

```
我要分析 R 內建的 AirPassengers（1949–1960 每月國際航空旅客數）。
這份要寫進統計諮詢報告，需要嚴謹、能交代統計依據，請幫我做時間序列分析。
```

Skill 會：

1. **反問你的分析取向**（模型的可解釋性 vs 預測導向）
2. 依你的選擇走對應路線，完成 EDA 到建模銜接
3. 產出 `report.qmd` + `report.html` + 乾淨資料

其他例子：

```
用 quantmod 抓台積電近三年股價，幫我做時間序列的初步分析。
我有每小時的伺服器流量，之後會持續進資料，想建可自動重跑的預測流程。
```

---

## 兩條路線的差別

| 面向 | (A) 模型的可解釋性 | (B) 預測準確 / pipeline 穩健 |
|---|---|---|
| 資料結構 | `ts` / `zoo` | `tsibble`（index + key）|
| 平穩性檢定 | **ADF**（不使用 KPSS）| **KPSS** + `unitroot_ndiffs()` 自動判定 |
| 缺失處理 | 明確的統計插補並說明理由 | `fill_gaps()` 先顯性化 |
| 定階 | ACF / PACF / EACF 人工判讀 | `ARIMA()` 自動定階 |
| 選模準則 | **BIC**（偏精簡、利於推論）| **AIC**（偏預測表現）|
| 繪圖 | base R | `ggplot2` / `feasts` |

> 兩種平穩性檢定的虛無假設**方向相反**（ADF 的 $H_0$ 是非平穩、KPSS 的 $H_0$ 是平穩），
> 報告中若同時引用務必註明，避免 p 值誤讀。

---

## 產出

| 檔案 | 說明 |
|---|---|
| `report.qmd` | Quarto 原始檔，可複核、修改、重跑 |
| `report.html` | 單一自足檔（CSS/JS/圖片內嵌），可直接寄送或列印 |
| 乾淨資料 | 已處理缺失／轉換，附前處理紀錄 |

報告格式規範：三線表（`flextable` + `theme_booktabs`）、中文標楷體 + 英數 Times New Roman、
圖說在下／表說在上、預測圖必附圖例、**預設不含 `sessionInfo()` 附錄**。

**預設不顯示程式碼**（程式碼仍在 `report.qmd` 內、編譯時照常執行）。需要含程式碼的版本時，跟 Claude 說「請產生含程式碼的報告」，或自行執行：

```
quarto render report.qmd -M echo:true -o report-with-code.html
```

---

## 檔案結構

```
timeseries-eda/
├── SKILL.md                          # 主流程（開場提問、五階段 EDA、報告規範）
├── assets/
│   └── report-template.qmd           # Quarto 報告範本（排版已設定）
├── references/
│   ├── classical-track.md            # (A) 傳統統計路線（含 ARIMAX、季節單根檢定）
│   ├── modern-fpp3-track.md          # (B) 現代 fpp3 路線（含動態迴歸）
│   ├── arima-garch-bridge.md         # 建模銜接：定階、auto.arima、ARCH/GARCH、Prophet
│   └── ljung-box-lags.md             # Ljung–Box 滯後期數 m 的依據（Tsay 講義查證）
└── evals/
    └── evals.json                    # 測試情境
```

採**漸進式載入**：`SKILL.md` 為精簡核心，`references/` 只在需要時才讀取。

---

## 參考資料

- Hyndman, R.J. & Athanasopoulos, G. (2021). *Forecasting: Principles and Practice*, 3rd ed. <https://otexts.com/fpp3/>
- Tsay, R.S. (2010). *Analysis of Financial Time Series*, 3rd ed.（Ch.2 迴歸模型含時間序列誤差；Ch.3 條件異質變異模型）
- tidyverts（tsibble / feasts / fable）<https://tidyverts.org>
- Quarto <https://quarto.org>
