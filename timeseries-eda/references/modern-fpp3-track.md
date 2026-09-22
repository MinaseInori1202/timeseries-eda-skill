# 現代 fpp3 路線（Modern track）

適用取向：**重視預測準確度與 pipeline 穩健／可重現**。
本路線以 Rob Hyndman《Forecasting: Principles and Practice (3rd ed.)》的 tidyverts 生態系為準。
出處：<https://otexts.com/fpp3/>。

## 生態系與套件

一次載入 `fpp3` 就會帶入下列核心套件（全部走 tidyverse 管線 `|>`）：

- **tsibble**：時間序列資料結構（`index` + `key` + 規則間隔）。
- **feasts**：視覺化與特徵萃取（`autoplot`、`gg_season`、`gg_subseries`、`gg_lag`、`ACF`、`STL`、`features`）。
- **fable**：建模與預測（本 skill 的 EDA 階段暫不深入，但差分階數自動判定會用到 `features()`）。
- 另含 `dplyr`、`ggplot2`、`lubridate` 等。

```r
library(fpp3)
```

## 1. 建立 tsibble

tsibble 的三個要素：**index（時間索引）**、**key（分組鍵，可多個）**、**interval（規則間隔）**。
即使沒 `select()`，tsibble 也會自動保留 index 與 key，以維持合法的時間結構。

依頻率選對時間類別函數：

| 頻率 | 函數 |
|------|------|
| 年 | 直接用整數年份 |
| 季 | `yearquarter()` |
| 月 | `yearmonth()` |
| 週 | `yearweek()` |
| 日 | `as_date()` / `lubridate::ymd()` |
| 次日（時/分/秒） | `as_datetime()` / `lubridate::ymd_hms()` |

```r
# 由既有 data frame 轉換
tsbl <- raw_data |>
  mutate(month = yearmonth(date)) |>
  as_tsibble(index = month, key = c(store, product))
```

### 隱含缺失值：先顯性化，不要隨手補

tsibble 的一大好處是能把「該有卻不見」的時間點抓出來。這些是 tsibble 套件本身的函數
（在 fpp3 正文各章不一定逐一示範，但屬 tidyverts 標準做法）：

```r
tsbl |> has_gaps()          # 有沒有間隔缺失
tsbl |> scan_gaps()         # 列出缺哪些時間點
tsbl |> count_gaps()        # 每條序列缺幾個
tsbl <- tsbl |> fill_gaps() # 把隱含缺失補成顯性的 NA（值先留 NA，之後再決定怎麼補）
```

`fill_gaps()` 只是把「跳號」補成一列值為 `NA` 的顯性缺失，**不是幫你填數值**。
要不要填、怎麼填（如 `tidyr::fill()` 前向填補、或 `imputeTS` 等）留到分解前再依資料性質決定，並在報告中記錄。

### 股價資料：先轉對數報酬

若是股價／價格型資產，依 SKILL.md 的通則，在建好 tsibble 後**先轉對數報酬再往下做**：

```r
tsbl <- tsbl |> mutate(log_return = difference(log(close)))
```

之後 `autoplot`、`ACF`、STL、KPSS 等都改用 `log_return`。可再畫 `tsbl |> ACF(log_return^2) |> autoplot()`
觀察波動聚集。

## 2. 視覺化探索

```r
tsbl |> autoplot(value)                       # 時間圖
tsbl |> gg_season(value)                      # 季節圖
tsbl |> gg_subseries(value)                   # 子序列圖
tsbl |> gg_lag(value, geom = "point")         # 落後圖
tsbl |> ACF(value) |> autoplot()              # 自相關圖（correlogram）
```

ACF 判讀：藍色虛線是顯著界線；有趨勢時小 lag 大且正、緩慢衰減；有季節時在季節倍數 lag 出現高峰；
趨勢+季節會呈現緩降且帶「扇貝狀」的形狀。

## 3. STL 分解

STL（Seasonal-Trend decomposition using Loess）比 X-11 / SEATS 泛用：可處理任何季節週期、
季節型態可隨時間演化、`robust = TRUE` 對離群值穩健。只做加法分解，乘法關係請先取對數。

```r
tsbl |>
  model(stl = STL(value ~ trend(window = 21) + season(window = "periodic"),
                  robust = TRUE)) |>
  components() |>
  autoplot()
```

- `trend(window=)`：控制趨勢平滑度（越大越平滑）。
- `season(window=)`：控制季節變動快慢；`"periodic"` 代表各週期季節型態固定。

## 4. 平穩性與差分（自動化）

現代路線用 **KPSS**（`unitroot_kpss()`，**H0 = 平穩**）搭配 **自動決定差分階數**，
而非人工反覆看圖猜：

```r
tsbl |> features(value, unitroot_kpss)     # KPSS：p < 0.05 → 需要差分
tsbl |> features(value, unitroot_ndiffs)   # 建議做幾次一階差分
tsbl |> features(value, unitroot_nsdiffs)  # 建議做幾次季節差分

# 差分（先季節差分，再一階差分）
tsbl <- tsbl |>
  mutate(d_value = difference(difference(value, lag = 12), lag = 1))
```

判讀與注意：

- KPSS 的 `p` 值被截在 0.01–0.1 之間；`p < 0.05` 代表**拒絕平穩**、需要差分。
- **順序**：先季節差分再一階差分——季節差分後有時就已平穩，不必再做一階。
- **過度差分**會誘發實際不存在的假自相關，`unitroot_ndiffs()` 的建議就是用來避免差過頭。

## 與傳統路線的對照（寫進報告時提醒讀者）

| 面向 | 現代 fpp3 | 傳統統計 |
|------|-----------|----------|
| 資料結構 | tsibble（index+key，時間語意由結構保證） | ts / zoo / xts 或排序好的 data.frame |
| 缺失處理 | `fill_gaps()` 先顯性化，再視需要填 | 直接統計插補（線性/季節/前向） |
| 平穩性檢定 | KPSS（H0=平穩） | ADF（H0=非平穩），兩者 p 值解讀相反 |
| 差分階數 | `unitroot_ndiffs()`/`nsdiffs()` 自動決定 | 人工看 ACF/PACF 與檢定反覆判斷 |
| 分解 | STL（穩健、季節可演化） | 傳統移動平均分解 / X-11 |
