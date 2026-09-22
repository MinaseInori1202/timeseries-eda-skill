# 傳統統計路線（Classical track）

適用取向：**重視可解釋性與統計嚴謹度**（寫論文、統計諮詢、需要為每個決策交代統計依據）。
本路線走經典的 Box–Jenkins / ARIMA 思維，以 **ADF 單位根檢定**與明確的統計轉換為主，
並在每一步都附上統計解讀與領域直覺（見 SKILL.md 的「每個結果都要解釋意義」原則）。

> 平穩性檢定本路線**只用 ADF，不使用 KPSS**。

## 生態系與套件

| 用途 | 套件 / 函數 |
|------|-------------|
| 時間序列物件 | `stats::ts()`、`zoo`、`xts` |
| 單位根檢定 | `tseries::adf.test()`（主）、`fUnitRoots::adfTest()` / `unitrootTest()`（可指定 type、較彈性）、`tseries::pp.test()`（Phillips–Perron，可交叉驗證）、`urca::ur.df()` |
| 缺失值插補 | `zoo::na.approx()`（線性）、`zoo::na.locf()`（前向填補）、`imputeTS`（季節性插補） |
| 變異數穩定 | `log()`（指數成長序列）、`forecast::BoxCox()` |
| 分解 | `stats::stl()`、`stats::decompose()` |
| 結構識別 | `stats::acf()`、`stats::pacf()`、`TSA::eacf()`（ARMA 定階）、`stats::ar()`（AIC 自動選階） |
| 差分 | `base::diff()`、`forecast::ndiffs()`、`forecast::nsdiffs()` |
| ARCH 效果檢定（金融序列） | `FinTS::ArchTest()`（Engle's LM test，檢定是否有波動聚集） |
| 殘差診斷（銜接建模） | `stats::Box.test(type = "Ljung-Box")`、`qqnorm()`/`qqline()`、`stats::shapiro.test()` |
| 波動度建模（銜接建模，超出本 skill） | `fGarch::garchFit()`（ARCH/GARCH 家族） |

本路線多用 R base graphics（`plot()`、`acf()`、`pacf()`），與統計諮詢報告的慣例一致；
若偏好 ggplot 風格，`forecast::ggAcf()` / `ggPacf()` 可替代。

```r
library(tseries)
library(zoo)
library(TSA)      # eacf
```

## 1. 建立 ts 物件

```r
# frequency 依季節週期設定：月 = 12、季 = 4、週 = 52、日(週季節) = 7
y <- ts(raw_data$value, start = c(2015, 1), frequency = 12)
```

若時間有跳號，先在 data.frame 階段把完整時間軸建好、對齊，再轉 `ts`——
**不要靠刪列來讓長度湊齊**，那會破壞時間間隔。

## 2. 時間圖與變異數穩定（log / Box–Cox）

先畫原始序列的時間圖，看趨勢、季節、循環與變異數是否隨水準放大：

```r
plot(y, col = "steelblue", lwd = 1.5,
     main = "Time Plot", ylab = "Value", xlab = "Year")
abline(h = mean(y), lty = 2, col = "gray50")
```

**判斷是否需要 log 轉換**：若序列呈**指數成長**、或季節波動幅度隨水準等比放大（乘法型季節），
先取對數把它轉成近似線性趨勢、加法型季節，穩定變異數後再往下做：

```r
log_y <- log(y)     # 指數成長 / 乘法季節時使用
```

> 解讀提醒：log 轉換不是為了好看，而是讓「變異數固定」這個平穩性前提更接近成立；
> 之後的差分、ADF、ARIMA 都在 log 尺度上進行，預測值最後再 `exp()` 轉回原尺度並說明。

（股價／價格型資產請改用對數報酬 `diff(log(price))`，見 SKILL.md 特例；
報酬序列可再用 `FinTS::ArchTest()` 檢定 ARCH 效果／波動聚集。）

## 3. 缺失值與異常值處理（明確的統計插補）

傳統路線在此就把缺失補好，並選一個「說得出統計理由」的方法：

```r
y <- na.approx(y)                 # 線性插值：適合平滑、無強季節的序列
# y <- na.locf(y)                 # 前向填補：適合階梯狀、狀態延續的序列
# y <- imputeTS::na_seasplit(y)   # 季節性插補：季節強時較合理
```

異常值先用時間圖與 STL 殘差找出來，是否處理要有領域理由，並在報告記錄。

## 4. 視覺化與結構識別（ACF / PACF / EACF / ar）

這是傳統路線判斷序列結構、為後續定階鋪路的核心。四樣工具搭配使用：

```r
par(mfrow = c(1, 2))
acf(as.numeric(y),  main = "ACF")
pacf(as.numeric(y), main = "PACF")
par(mfrow = c(1, 1))

TSA::eacf(as.numeric(y))          # 擴展 ACF，找 ARMA(p,q) 候選頂點
ar(y, method = "mle")             # 以 AIC 自動選 AR 階數（可解釋性參考）
```

**判讀規則（Box–Jenkins）**：

- ACF 拖尾、PACF 在 p 階後截斷 → 偏 **AR(p)**。
- ACF 在 q 階後截斷、PACF 拖尾 → 偏 **MA(q)**。
- 兩者都拖尾 → 偏 **ARMA**，用 **EACF** 找「o 三角形頂點」定 (p, q)，再以 AIC/BIC 裁決。
- 有季節時，在季節 lag（12、24… 或 4、8…）觀察季節性的 AR/MA 徵兆。

`ar()` 常選出較高階數（反映序列持續性）；報告時可**兼顧可解釋性**，挑選有明確圖形依據的較精簡階數，並說明取捨理由。

（本 skill 以 EDA 為範圍；上面的識別結果是要交給後續 ARIMA 建模用，不在此配適模型。）

## 5. 分解

```r
decomp <- stl(y, s.window = "periodic", robust = TRUE)
autoplot(decomp)                  # 或 plot(decomp)
```

`s.window = "periodic"` 表季節型態固定；需要季節隨時間演化時改成具體的奇數視窗長度。
逐一解讀趨勢、季節、殘差三成分：殘差是否仍有結構？有的話代表尚有未被解釋的規律。

## 6. 平穩性與差分（以 ADF 為中心）

以 **ADF（Augmented Dickey–Fuller）** 判定平穩，**H0 = 有單位根（非平穩）**：

```r
adf.test(y)          # 或 adf.test(log_y)；p < 0.05 → 拒絕 H0 → 判定平穩
# pp.test(y)         # Phillips–Perron，另一種單位根檢定，可交叉驗證
```

- **解讀方向**：ADF 的 `p < 0.05` 代表**平穩**。
- 嚴謹一點可用 `urca::ur.df()` 指定 `type = c("none", "drift", "trend")`，看 τ 統計量與臨界值，
  以區分「含漂移／趨勢」的情形，並在報告說明選了哪一型與理由。
- 不平穩就差分，用 `ndiffs()` / `nsdiffs()` 輔助判斷次數：

```r
ndiffs(y)                                    # 建議一階差分次數
nsdiffs(y)                                   # 建議季節差分次數
y_d <- diff(diff(log_y, lag = 12), differences = 1)   # 先季節差分，再一階差分
adf.test(y_d)                                # 差分後再檢定，確認已平穩
```

- **順序**：先季節差分（`lag = 頻率`，季資料 4、月資料 12）再一階差分。
- **過度差分**要避免——差過頭會製造假自相關、放大變異數。差分後重畫 ACF/PACF 確認結構，並決定 (d, D)。

## 銜接後續建模（提示，不在本 skill 範圍）

EDA 完成後，傳統路線通常接著配適 ARIMA/SARIMA、比較 AIC/BIC、並做**殘差診斷**才算完整：

- 殘差白化檢定：`Box.test(resid, lag = 20, type = "Ljung-Box", fitdf = 參數個數)`，p > 0.05 表殘差已白化。
- 殘差常態性：`qqnorm()` + `qqline()`、`shapiro.test()`。
- 預測附 95% 預測區間，log 尺度預測記得 `exp()` 轉回並解讀。

先把上述識別與平穩化結果整理清楚，讓建模階段能直接接手。

## 產出提醒

本路線的賣點是「可解釋」，報告裡每個決策都要寫清楚**用了什麼檢定、統計量與 p 值、為何這樣轉換／補值、
為何差分幾次**，並附上領域直覺解讀，讓讀者（或審稿人）能複核。
