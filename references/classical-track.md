# 傳統統計路線（Classical track）

適用取向：**重視模型的可解釋性**（寫論文、統計諮詢、需要為每個決策交代統計依據）。
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
| 有外生變數（ARIMAX，第 7 節） | `forecast::Arima(xreg=)` / `auto.arima(xreg=)`、`TSA::arimax()`、`lmtest`（`dwtest`/`bgtest`/`coeftest`）、`tseries::po.test()`、`urca::ca.jo()`、`sandwich::NeweyWest()` |
| 殘差診斷（銜接建模） | `stats::Box.test(type = "Ljung-Box")`、`qqnorm()`/`qqline()`、`stats::shapiro.test()` |
| 波動度建模（銜接建模，超出本 skill） | `rugarch`（GARCH/EGARCH，金融常用）、`fGarch::garchFit()`（ARCH/GARCH 家族） |

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

> **長季節／日資料的年週期要特別處理（重要）。** 日資料常是**年週期 s≈365**（甚至同時有週季節 7），
> 這種**長季節**或**多重季節****不要用 `ts(frequency = 365)` 去配 SARIMA(s=365)**——階數太高、幾乎估不動、也不穩。
> 改用下列任一：
> - **STL 分解**先看季節（`stl` / `mstl`），或用 `forecast::stlf()` 直接做「去季節後建模再貼回」的預測；
> - **動態諧波迴歸（dynamic harmonic regression）**：把年季節用少數傅立葉諧波項顯式建模、誤差走非季節 ARIMA：
>   ```r
>   library(forecast)
>   y  <- msts(raw_data$value, seasonal.periods = c(7, 365.25))   # 多重季節
>   xr <- fourier(y, K = c(3, 6))                                  # 每個週期的諧波階數 K
>   fit <- auto.arima(y, xreg = xr, seasonal = FALSE)             # 季節交給諧波、ARIMA 只管殘差
>   fc  <- forecast(fit, xreg = fourier(y, K = c(3, 6), h = 30))
>   ```
>   `K` 由小到大試、以 AICc 最小者為準（K 越大季節曲線越彎）。這也是現代路線 `fable` 用 `fourier()` 的同一思路。

## 2. 時間圖與變異數穩定（log / Box–Cox）

先畫原始序列的時間圖，看趨勢、季節、循環與變異數是否隨水準放大：

```r
plot(y, col = "steelblue", lwd = 1.5,
     main = "", ylab = "數值（單位）", xlab = "時間")   # 不放圖內主標題，圖名由 fig-cap 在下方給
```

**判斷是否需要 log 轉換**：若序列呈**指數成長**、或季節波動幅度隨水準等比放大（乘法型季節），
先取對數把它轉成近似線性趨勢、加法型季節，穩定變異數後再往下做：

```r
log_y <- log(y)     # 指數成長 / 乘法季節時使用
```

**判斷時要用圖與量化數字（不要只憑感覺）**：先看原始尺度的時間圖（上面的 `plot(y)`）與季節圖（`forecast::ggseasonplot(y, continuous = TRUE)`），再量化逐年年內全距相對於年平均的放大倍數，轉換前後各算一次：

```r
yr <- floor(time(y))                                   # 月、季資料的年份
rng <- tapply(y, yr, function(z) diff(range(z)))       # 逐年年內全距
mu  <- tapply(y, yr, mean)                             # 逐年平均水準
rbind(range = rng, mean = mu)                          # 全距放大速度 > 水準放大速度 → 乘法型、需要轉換
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
# 最方便：一次畫出「序列圖 + ACF + PACF」，定階時很好用
forecast::ggtsdisplay(y, main = "")

# 或分開畫：
par(mfrow = c(1, 2))
acf(as.numeric(y),  main = "自相關函數 ACF")
pacf(as.numeric(y), main = "偏自相關函數 PACF")
par(mfrow = c(1, 1))

TSA::eacf(as.numeric(y))          # 擴展 ACF，找 ARMA(p,q) 候選頂點
ar(y, method = "mle")             # 以 AIC 自動選 AR 階數（可解釋性參考）
```

**判讀規則（Box–Jenkins）**：

- ACF 拖尾、PACF 在 p 階後截斷 → 偏 **AR(p)**。
- ACF 在 q 階後截斷、PACF 拖尾 → 偏 **MA(q)**。
- 兩者都拖尾 → 偏 **ARMA**，用 **EACF** 找「o 三角形頂點」定 (p, q)，再以 AIC/BIC 裁決。
- 有季節時：**先依週期長度做季節差分再定階**，並分「季節部分 (P,Q)：只看 s,2s,3s 週期點」與
  「非季節部分 (p,q)：只看前半段短期 lag」兩段判讀——完整規則見 SKILL.md 的「SARIMA 定階規則」。

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

### 季節單根檢定（有季節性時，決定 D 的依據）

標準 ADF 只檢定**零頻率（非季節）單根**，對季節頻率的單根不敏感，也可能被確定性季節誤導（報告中常見：未差分序列 ADF 也拒絕 H₀）。
因此月、季資料要用**季節單根檢定**決定季節差分階數 D，**H₀ 都是「存在季節單根（需要季節差分）」**，p 小或統計量低於臨界值 → 拒絕、不需季節差分：

```r
library(forecast)
nsdiffs(y, test = "ocsb")      # OCSB（Osborn–Chui–Smith–Birchenhall）：回傳建議的 D
forecast:::ocsb.test(y)        # 取得統計量與 5% 臨界值（內部函數，供報告引用數字）
# HEGY（逐季節頻率檢定）：需安裝 uroot
# uroot::hegy.test(y, deterministic = c(1, 1, 1), lag.method = "AIC", maxlag = 12)
```

- **不要用預設的 `nsdiffs(y)` 當「檢定」**：預設 `test = "seas"` 是 STL 季節強度閾值（強度 > 0.64 就建議差分），不是單根檢定。
- **`test = "ch"`（Canova–Hansen）H₀ 是「季節型態平穩」，方向與 KPSS 相同，本路線不使用**（原則同「不用 KPSS」）。
- **HEGY**：H₀ 為各季節頻率單根；月資料會輸出多個頻率的統計量（零頻率、π、各季節諧波…），逐一判讀，只有部分頻率有單根時要說明（可能只需部分季節差分或季節虛擬變數）。**`uroot` 套件尚未實測，使用前須安裝。**
- **結論不一致時（實測：log(AirPassengers)，季節強度建議 D = 1，OCSB 統計量 −1.95 < 5% 臨界值 −1.80 → 拒絕季節單根 → D = 0；原始尺度則兩者皆 D = 1）**：
  兩個結論都要寫進報告，說明取捨依據（時間圖／季節圖、差分後 ACF 在季節落後期的表現、過度差分風險、殘差診斷），**不可只挑與預期相符的那個**。最終取捨由使用者依分析目的決定。
- 報告中的檢定結果表要列出 H₀、統計量、臨界值或 p 值、判定。
- **過度差分**要避免——差過頭會製造假自相關、放大變異數。差分後重畫 ACF/PACF 確認結構，並決定 (d, D)。

## 7. 有外生變數：迴歸模型含時間序列誤差（ARIMAX）

**何時用**：除了 y 本身，還有其他可觀測的解釋變數 x（如氣溫、利率、廣告量、1 月虛擬變數…），
想估計 x 對 y 的影響，或用 x 提升預測。依據 Tsay, *Analysis of Financial Time Series* Ch.2「Regression Models with Time Series Errors」。

**模型定義（務必寫進報告）**：採「**迴歸誤差項為 ARMA**」的參數化——

> y<sub>t</sub> = β₀ + β₁x<sub>t</sub> + W<sub>t</sub>，W<sub>t</sub> 服從 ARMA(p, q)。

這正是 `forecast::Arima(xreg = )`、`stats::arima(xreg = )`、`fable::ARIMA(y ~ x)` 的做法，β 的詮釋為「x 每增加 1 單位，y 的水準平均變動 β」。
**不是** y<sub>t</sub> = β₀ + β₁x<sub>t</sub> + φ₁y<sub>t−1</sub> + …（把 y 的落後項直接放進迴歸）：後者的 y<sub>t−1</sub> 已含 x<sub>t−1</sub> 的效應，β 的意義會被混淆。

**為什麼不能直接 `lm()` 就解釋**：誤差若有序列相關，OLS 的標準誤與 t 值不可靠（常被低估而灌水）；
兩條都不平穩的序列互相迴歸，還會產生**假性迴歸（spurious regression）**——R² 很高、t 很顯著，卻毫無因果或長期關係。

### 三步驟流程與每一步的 R 檢查

| Tsay 步驟 | 要做的事 | R 套件／函數 |
|---|---|---|
| **① 先配適線性迴歸，檢查殘差序列相關** | `lm()` 後檢查殘差 | `lmtest::dwtest()`（Durbin–Watson，僅 lag-1）、`lmtest::bgtest(order = k)`（Breusch–Godfrey，任意階）、`Box.test(resid, type = "Ljung-Box")`（**講義建議以 Ljung–Box 取代 DW**）、`forecast::ggtsdisplay(resid)` |
| **② 殘差是否單根非平穩？** | **是** → y 與 x **都**取一階差分，回到 ① 重做；**否** → 替殘差定 ARMA 階數 | 單根：`tseries::adf.test(resid)`（僅作參考，見下方警告）、長期關係用 `tseries::po.test(cbind(y, x))`（Phillips–Ouliaris 共整合檢定）或 `urca::ca.jo()`（Johansen）；定階：`ggtsdisplay(resid)`、`TSA::eacf(resid)`、`auto.arima(resid, d = 0, seasonal = FALSE, ic = "bic")` |
| **③ 最大概似聯合估計並檢查** | 迴歸係數與 ARMA 參數**一起**估計，再診斷、改進 | `forecast::Arima(y, order = c(p,d,q), xreg = x)` 或 `forecast::auto.arima(y, xreg = x)`；`TSA::arimax()`（亦支援轉換函數）；係數顯著性 `lmtest::coeftest(fit)`；殘差 `checkresiduals(fit)`、`Box.test(resid(fit), lag, type = "Ljung-Box", fitdf = p + q)` |

```r
library(forecast); library(lmtest); library(tseries); library(TSA)

## ① 先配適線性迴歸並檢查殘差
m0 <- lm(y ~ x)                                  # y、x 為 ts 物件
dwtest(m0)                                       # DW ≈ 2(1 − ρ̂₁)；遠小於 2 → 正自相關
bgtest(m0, order = 8)                            # 較高階的序列相關檢定（季資料 8、月資料 12 起試）
Box.test(resid(m0), lag = 10, type = "Ljung-Box")
ggtsdisplay(resid(m0), main = "")

## ② 殘差是否仍有單根？（決定要不要差分 y 與 x）
po.test(cbind(y, x))                             # H0：無共整合（殘差有單根）；p < 0.05 才能說 y 與 x 有長期均衡關係
## 殘差非平穩且無共整合（po.test p 大）→ 假性迴歸，y 與 x 一起差分，回到 ①
## 有共整合但殘差序列相關仍強（如下方利率例）→ 講義仍選差分，聚焦「變動對變動」的短期關係；
##   此時 β 詮釋改為變動的效應、並放棄水準層次的長期資訊，須在報告說明
cy <- diff(y); cx <- diff(x)
m1 <- lm(cy ~ cx - 1)                            # 差分後的迴歸（漂移項視需要保留）
Box.test(resid(m1), lag = 10, type = "Ljung-Box")
## 若殘差平穩 → 替殘差定 ARMA 階數
TSA::eacf(resid(m0))
auto.arima(resid(m0), d = 0, seasonal = FALSE, ic = "bic")   # ic 依開場取向：解釋 BIC／預測 AIC

## ③ 聯合最大概似估計 + 檢查
fit <- Arima(y, order = c(1, 0, 2), xreg = x)    # 階數由 ② 決定；或直接 auto.arima(y, xreg = x, seasonal = FALSE, ic = "bic")
coeftest(fit)                                    # 迴歸係數與 ARMA 參數的 z 檢定
checkresiduals(fit)                              # 殘差圖 + Ljung-Box（自動扣掉模型自由度）
Box.test(resid(fit), lag = 10, type = "Ljung-Box", fitdf = 3)   # fitdf = AR + MA 參數個數（不含迴歸係數）

## 解釋力：差分、無截距（include.mean = FALSE）的模型用未中心化平方和算 R²（Tsay 講義第 60 頁）
rsq <- (sum(cy^2) - sum(resid(fit)^2)) / sum(cy^2)   # cy 為差分後的 y；只適用無截距的差分模型
```

> **警告 1｜OLS 殘差的 ADF 檢定不可靠。** 殘差是估計出來的，標準 ADF 臨界值不適用（會過度拒絕，
> 把「假性迴歸」誤判成平穩）。實測：兩條獨立隨機漫步的 OLS 殘差，`adf.test()` p = 0.03，
> 但 `po.test()` p = 0.15——後者才正確。判斷是否需要差分，**以 `po.test()` / Johansen 檢定為準**，
> `adf.test(resid)` 只當輔助、且須在報告註明這個限制。
>
> **警告 2｜只用 OLS 殘差定階會失準。** 步驟 ② 由 OLS 殘差得到的階數只是**起點**；
> 實測 OLS 殘差的 BIC 選了白噪音，聯合估計後卻是 ARMA(1,2) 且各參數顯著，β 也從 0.27 變成 0.20。
> 所以步驟 ③ 的聯合估計與殘差檢查**不能省**，必要時直接用 `auto.arima(y, xreg = x)` 與 ② 的人工階數比較。
>
> **警告 3｜序列相關不能只看 DW。** DW 只看 lag-1；報告以 **Ljung–Box／Breusch–Godfrey** 為主，DW 只當補充。

> **警告 4｜`arima()` 與 `Arima()` 的 AIC 不要混用比較。** 同一模型（同樣的概似值 3136.62），R 4.5.1 的 `stats::arima()` 印出 AIC = −6269.23，
> `forecast::Arima()` 則是 −6267.23（講義印的也是後者）。比較模型時**同一批用同一個函數**。

### 實例驗證：Tsay 講義的利率資料（w-gs1yr.txt／w-gs3yr.txt，2467 週）

講義第 58–60 頁的數字全部重現（`lm` 斜率 0.92955、R² = 0.9649；差分後斜率 0.791935；聯合估計 MA(1) θ = 0.1823、β = 0.7936、概似值 3136.62；R² = 0.8310），且 `Arima(r3, order = c(0,1,1), xreg = r1)` 在水準值上直接處理與手動差分結果完全相同。
照本節流程多做的檢查，則得到講義沒有的結論：

| 檢查 | 結果 | 意義 |
|---|---|---|
| ① 水準迴歸殘差 | DW = 0.020、Ljung–Box／BG p < 2e-16 | 嚴重序列相關，OLS 的 t 值不可信 |
| ② `po.test` | p < 0.01（有共整合） | 利率間有長期均衡關係，**不是**假性迴歸；講義仍選差分以處理短期變動 |
| ③ 差分後 `Box.test(resid(fit), 12, fitdf = 1)` | p ≈ 1e-7 | **MA(1) 並未榨乾序列相關**（講義只看 ACF 目視）；`auto.arima(ic = "aic")` 傾向 ARMA(5,4)，`"bic"` 則回到 MA(1) |
| 第二層：`FinTS::ArchTest(resid(fit), lags = 12)` | p < 2e-16 | 強烈 ARCH 效應，標準誤與區間需考慮波動聚集（接 GARCH） |
| HAC | x 的標準誤 0.0073 → 0.0203 | 序列相關／異質變異下，原始 t 值被高估約 2.8 倍 |

報告中應如實寫出：MA(1) 是 BIC 下的精簡模型，但 Ljung–Box 與 ARCH 檢定顯示尚有殘留結構，結論與預測區間須保守解讀。

### 差分、截距與預測的細節

- `Arima(y, order = c(p,1,q), xreg = x)`：**d ≥ 1 時 xreg 會連同 y 一起差分**，截距項自動消失，係數是「x 的變動對 y 的變動」的效應。
- `include.mean = FALSE`（講義例）：y 已中心化或 x 含虛擬變數時可不要截距。
- **預測需要未來的 x**：`forecast(fit, xreg = x_future)`，`x_future` 列數 = h、欄位順序同配適時。
  x_future 可來自 x 自己的預測、已知排程（假日、廣告檔期）或情境假設（如「維持歷史平均」），
  **務必在報告說明來源**；用 x 的預測值時，預測區間**不含 x 本身的不確定性**，要另行註明。
- **虛擬變數／季節 x**：月份、假日、介入事件可做成 0/1 欄位放進 `xreg`（講義第 52–53 頁用 `jan` 虛擬變數的例子）。
- **殘差平方檢定**：若 ARIMAX 殘差仍有 ARCH（金融資料常見），均值方程用 `xreg`、波動度接 GARCH——
  `rugarch::ugarchspec(mean.model = list(armaOrder = c(p, q), external.regressors = as.matrix(x)))`，見 `arima-garch-bridge.md`。
- **只要結論不依賴 ARMA 誤差（例如只做係數推論）**：可用 OLS + 穩健標準誤 `lmtest::coeftest(m0, vcov. = sandwich::NeweyWest(m0, lag = 8, prewhite = FALSE))`（HAC，Newey–West）；
  實測 us_change 資料中 x 的標準誤從 0.047 變為 0.076，t 值顯著下降。注意這只修正標準誤，**不改善預測**。

### 報告寫法

交代：(1) 使用的模型定義（迴歸誤差項為 ARMA）；(2) ① 的序列相關檢定與 p 值；(3) ② 是否差分及依據（含警告 1 的限制）；
(4) ③ 的最終模型、β 的估計值與標準誤／z 檢定、對 x 的**白話詮釋**；(5) 殘差兩層診斷（見 SKILL.md）；(6) 預測所採用的 x 未來值假設。

## 銜接後續建模（提示，不在本 skill 範圍）

EDA 完成後，傳統路線通常接著配適 ARIMA/SARIMA、依準則比較、並做**殘差診斷**才算完整：

- **選模準則先問使用者**：預測導向用 **AIC**、解釋／統計推論導向用 **BIC**（BIC 懲罰較重、傾向較精簡的模型）。
  傳統路線多屬 (A) 可解釋取向，常以 BIC 為主，但仍應向使用者確認。


- 殘差白化檢定：`Box.test(resid, lag = 10, type = "Ljung-Box", fitdf = p + q)`，p > 0.05 表殘差已白化。
- 殘差常態性：`qqnorm()` + `qqline()`、`shapiro.test()`。

**預測：點估計 + 預測區間圖（務必畫出來）**。用 `forecast` 套件最省事——`forecast()` 產生點預測與 80%／95% 預測區間，`autoplot()` 直接把預測疊在歷史序列上：

```r
library(forecast); library(ggplot2)
fc <- forecast(fit, h = 12)               # fit 為 Arima()/auto.arima() 物件；h = 預測期數
autoplot(fc) +
  autolayer(fc$mean, series = "點預測") +          # 讓「點預測」進圖例
  labs(title = "預測：點估計與 80% / 95% 預測區間",
       x = "時間", y = "數值",
       colour = "序列", fill = "預測區間") +        # 區間層級(80/95)列在 fill 圖例、點預測列在 colour 圖例
  theme(legend.position = "right")
# 取數值：fc$mean（點預測）、fc$lower / fc$upper（各含 80%、95% 兩欄）
```

> **一定要有圖例**：`forecast::autoplot(fc)` 的陰影區間預設會有 `level`（80/95）圖例；再用
> `autolayer(fc$mean, series = "點預測")` 讓點預測也出現在圖例，讀者才分得清點估計與區間。

- **log／Box–Cox 尺度**：用 `Arima(y, lambda = 0)` 或 `forecast(fit, lambda = 0, biasadj = TRUE)`，`autoplot` 會**自動反轉換回原尺度**、區間也一併還原（不用自己 `exp()`）。
- **手動法**（如作業風格，base plot）：`p <- predict(fit, n.ahead = 12); upper <- p$pred + 1.96*p$se; lower <- p$pred - 1.96*p$se`，
  再把歷史末點接到預測首點畫線、加上區間虛線，並用 `legend()` 標出「觀測值 / 點預測 / 95% 區間」：
  ```r
  legend("topleft", legend = c("觀測值", "點預測", "95% 區間"),
         col = c("black", "tomato", "tomato"), lty = c(1, 1, 2), lwd = c(1.5, 2, 1), bty = "n")
  ```
- **判讀**：預測區間會**隨 horizon 變寬**（不確定性累積）；金融報酬的點預測常接近常數，重點在區間寬度與波動（條件變異數／GARCH 的區間更有意義）。

先把上述識別、平穩化、預測結果整理清楚，讓建模階段能直接接手。

## 產出提醒

本路線的賣點是「可解釋」，報告裡每個決策都要寫清楚**用了什麼檢定、統計量與 p 值、為何這樣轉換／補值、
為何差分幾次**，並附上領域直覺解讀（報告本文以客觀敘述呈現，不寫「以便讀者複核」之類字句）。
