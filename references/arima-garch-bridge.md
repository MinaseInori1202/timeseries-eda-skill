# 銜接建模：SARIMA 定階 → 自動選階 → 殘差診斷 → ARCH/GARCH

本檔是「EDA 完成、平穩性確認後」的**建模銜接**細節。實際配適模型超出本 skill 的 EDA 範圍，
但「該用什麼階數、要不要上 GARCH」這些**判斷與方向**要在報告中點明。
前置：已完成 SKILL.md 的第 5 階段（平穩性、差分、畫出 ACF/PACF 並做基本判讀）。

## SARIMA 定階規則（有週期波動時務必照這個順序）

- **先處理季節、再定階**：若序列有週期 `s` 的波動（月資料 s=12、季資料 s=4…），
  **先依週期長度做季節差分**（`diff(x, lag = s)`）讓序列平穩，**再**用下面規則定階；否則季節結構會干擾 ACF/PACF 判讀。
  - **例外：長季節／多重季節（日資料年週期 s≈365 等）不適用 SARIMA。** 季節差分 lag=365、SARIMA(s=365) 幾乎估不動也不穩。
    改用**動態諧波迴歸（Fourier 諧波項 + 非季節 ARIMA 誤差）**或 **STL**——見 classical-track 的「長季節」段與 modern-fpp3-track 的 `fourier()` 範例。
- **季節部分 `(P, Q)`——只看季節週期點 `s, 2s, 3s…`**：
  - 季節 MA `Q`：ACF 在 lag `Qs` 截斷（例如 lag `s` 顯著、`2s` 之後立刻落回陰影區），且 PACF 緩慢衰減 → 設 `Q = 1`。
  - 季節 AR `P`：PACF 在 lag `Ps` 截斷，且 ACF 緩慢衰減 → 設 `P = 1`。
  - 實務：金融與經濟資料的季節階數 `P`、`Q` 很少超過 1 或 2。
- **非季節部分 `(p, q)`——只看前半段短期 lag**：
  - MA `q`：ACF 在 lag `q` 截斷，PACF 拖尾。
  - AR `p`：PACF 在 lag `p` 截斷，ACF 拖尾。

## 全自動選項：auto.arima（可搭配手動判讀一起用）

除了手動看圖定階，最著名的自動函數是 `forecast::auto.arima()`，
它會**網格搜索、依 AIC/BIC 自動選出非季節 `(p,d,q)` 與季節 `(P,D,Q)` 的最佳階數**。
可把**已經決定好的差分次數鎖定**（如季節差分 `D = 1`、一般差分 `d = 1`），讓它只搜尋其他階數；
並用 `ic = "aic"` / `"bic"` 對應選模取向（預測→AIC、解釋/推論→BIC）：

```r
forecast::auto.arima(y, d = 1, D = 1, ic = "bic",
                     stepwise = FALSE, approximation = FALSE)  # 關閉近似做較完整搜尋
```

即使用自動選階，**仍建議先手動畫 ACF/PACF 判讀**——理解資料結構，並驗證自動選出的階數是否合理，不要全盤照收。

## 推薦的效率流程（先自動、殘差沒過再手動）

1. 用 `ts()` 設好週期（`frequency`），丟進 `auto.arima()` 跑出一個**基準模型**。
2. 執行 `forecast::checkresiduals(fit)`（會一併做**殘差的 Ljung-Box 檢定**）。
3. 依殘差 Ljung-Box 決策：
   - **p 值 > 0.05**（殘差已白化、沒有殘留自相關）→ **均值方程（mean equation）**已足夠，不必再看 EACF；
     **非金融資料**可到此收工。**但金融／報酬資料不要停手**——一階矩沒有自相關 ≠ 沒有結構，
     接著 **100% 要檢「殘差平方的 ARCH 效應」**（見下方「殘差平方 ARCH」）。
   - **p 值 < 0.05**（殘差還有沒榨乾的資訊）→ 這時才手動介入，用 `ACF/PACF` 或 `TSA::eacf()` 嘗試**微調（通常調高階數）**。
4. 所以 **`EACF` 不是必要步驟，而是「`auto.arima` 選出的模型沒通過殘差檢定」時的輔助工具**——
   用來觀察是否有被忽略的結構，再回頭調整階數。

> 注意區分兩種 Ljung-Box：股價分支是對**報酬序列本身**做（決定走 ARMA 或 ARCH）；
> 這裡是對**配適後的殘差**做（判斷模型是否足夠）。兩者用途不同，別混淆。

## 殘差平方的 ARCH 效應檢驗（金融資料必做）

**即使殘差通過一階 Ljung-Box（均值特徵已捕捉完），金融資料也不能就此結束。** 金融數據的核心本質是：
**報酬率／殘差本身（一階矩）可能已無自相關，但其「平方」（二階矩）通常仍有強烈自相關**（波動度聚集）。
`checkresiduals()` 只檢一階、不會自動檢殘差平方，所以要額外做：

1. **視覺化**：畫**殘差平方**的 `ACF` 與 **`PACF`**——`acf(residuals(fit)^2)`、`pacf(residuals(fit)^2)`。
   有 ARCH 效應時，平方後自相關訊號會全部浮現；**殘差平方的 PACF 截斷位置可初判 ARCH 階數 m**。
2. **正式檢定（Tsay 給兩種，可並用）**：
   - **Ljung-Box Q(m) 檢殘差平方**：`Box.test(residuals(fit)^2, lag = 12, type = "Ljung")`。
   - **Engle 的 LM（Lagrange Multiplier）檢定**：`FinTS::ArchTest(residuals(fit), lags = 12, demean = FALSE)`
     （注意 `ArchTest` 作用在殘差本身、不是平方；`demean` 控制要不要先去均值）。
   兩者 **p < 0.05** 都代表有顯著 ARCH 效應。
   - **p < 0.05** → 拒絕虛無假設，證實殘差有顯著 **ARCH 效應** → **升級到 GARCH**：
     保留原本的 ARIMA／SARIMA 當作 **Mean Equation**，在其上疊加 **GARCH(1,1)**（或 EGARCH）動態建模波動度；
     R 中常用 **`rugarch`** 套件（`fGarch` 亦可）。
   - **p ≥ 0.05** → 沒有顯著 ARCH，波動度可視為固定。

## GARCH 配適與終極診斷（rugarch）

1. 用 `rugarch::ugarchspec()` 定義 **ARMA/SARIMA + GARCH(1,1)** 結構：`mean.model` 設 `armaOrder`（均值方程，
   季節項通常另以外生變數或先前差分處理），`variance.model` 設 `garchOrder = c(1, 1)`。
2. 把**原始序列（不要自己先差分）丟進 `ugarchfit()`**——rugarch 會依規格**自動處理差分與均值**。
3. `show(fit)`，拉到輸出最下方看**兩大核心指標**確認模型是否完善：
   - **Optimal Parameters（參數顯著性）**：檢查 `ar1, ma1, alpha1, beta1` 等所有參數 p 值是否都 **< 0.05**；
     都顯著代表階數選得很健康。
   - **Robust Standard Errors 下的兩個加權 Ljung-Box（終極殘差檢驗）**：
     - **Weighted Ljung-Box on Standardized Residuals**：檢查「均值」有沒有榨乾，**p 必須 > 0.05**。
     - **Weighted Ljung-Box on Standardized Squared Residuals**：檢查「波動度（GARCH 效果）」有沒有榨乾，**p 也必須 > 0.05**。

```r
library(rugarch)
spec <- ugarchspec(
  mean.model         = list(armaOrder = c(p, q), include.mean = TRUE),
  variance.model     = list(model = "sGARCH", garchOrder = c(1, 1)),
  distribution.model = "std"          # 報酬厚尾，常用 t 分布
)
fit <- ugarchfit(spec, data = returns)  # 丟原始（報酬）序列，rugarch 自動處理均值/差分
show(fit)
```

**當兩個加權 Ljung-Box 的 p 值都 > 0.05**：代表不論是漲跌趨勢（均值）或風險波動（變異數），
都已被模型完整捕捉——分析才算真正完善。

## 波動度建模的標準程序與模型家族（Tsay《AFTS》Ch.3）

**建立波動度模型的五步驟**（與均值的 Box–Jenkins 對稱）：

1. **建均值方程並檢定 ARCH**：先配 ARMA／SARIMA 去掉線性相關，再對殘差檢 ARCH 效應（上面兩種檢定）。
2. **定階**：用**殘差平方的 PACF** 判斷 ARCH 階數 m；GARCH 則多從 **GARCH(1,1)** 起手（實務最常足夠）。
3. **估計**：條件最大概似（conditional MLE），均值與波動度**一起估**。
4. **模型檢查**：用**標準化殘差** `ã_t = a_t / σ_t`——`Ljung-Box(ã_t)` 檢均值是否足夠、`Ljung-Box(ã_t²)` 檢波動度是否足夠
   （都要 p > 0.05）；再用 `ã_t` 的**偏態、峰態、QQ 圖**檢查分布假設。（這正是 `rugarch::show(fit)` 印出的那幾項。）
5. **軟體**：R（`rugarch` 或 `fGarch`）。

**誤差分布（厚尾）**：報酬多厚尾，`ε_t` 別只用常態；可選 **Student-t、skewed-t、GED**。
rugarch：`distribution.model = "std"／"sstd"／"ged"`；fGarch：`cond.dist = "std"／"sstd"／"ged"`。

**模型家族（GARCH(1,1) 不夠時往這裡挑）**：

- **EGARCH／TGARCH(GJR)**：捕捉**槓桿效應／不對稱**——大跌比大漲更會推升波動（`σ_t²` 對正、負衝擊反應不同）。
  rugarch `model = "eGARCH"`／`"gjrGARCH"`。
- **IGARCH**：`α1 + β1 ≈ 1`、波動高度持續（衝擊幾乎不消退）。
- **GARCH-M**：把波動度放進均值方程（風險溢酬，報酬隨風險上升）。

**fGarch 寫法（Tsay 課堂用法，可替代 rugarch）**：

```r
library(fGarch)
m <- garchFit(~ arma(1, 0) + garch(1, 1), data = returns,
              cond.dist = "std", trace = FALSE)   # cond.dist：norm/std/sstd/ged
summary(m)     # 含係數顯著性、標準化殘差 Ljung-Box(R 與 R^2)、LM Arch、Jarque-Bera、資訊準則
```

**平穩與持續性**：GARCH(1,1) 弱平穩需 `α1 + β1 < 1`；越接近 1 越持續（趨向 IGARCH）。
**波動度預測**：多步預測 `σ_t²(h)` 會隨 h **收斂到無條件變異數** `α0 / (1 − α1 − β1)`。

### 自動產生模型方程式（LaTeX；取代 equatiomatic）

`equatiomatic` 不支援 ARIMA/GARCH。但 ARIMA/GARCH 的**方程式框架是固定的**
（均值 `r_t = μ + Σφ_i r_{t-i} + … + ε_t`；波動度 `σ_t² = ω + α ε_{t-1}² + β σ_{t-1}²`），
實務上寫一個小函式**讀取估計係數、用 `sprintf()` 填進 LaTeX 模板**即可，讓報告直接印出帶實際數字的方程式：

```r
extract_garch_eq <- function(fit) {
  coefs <- fit@fit$coef                    # rugarch fit 的估計係數
  # 依實際模型調整要抓的名稱（先看 names(coefs)）；此處以 AR(1)+GARCH(1,1) 為例
  mu    <- round(coefs["mu"],     4)
  ar1   <- round(coefs["ar1"],    4)
  omega <- round(coefs["omega"],  4)
  alpha <- round(coefs["alpha1"], 4)
  beta  <- round(coefs["beta1"],  4)

  sprintf("
\\begin{aligned}
r_t &= %s + %s\\, r_{t-1} + \\epsilon_t \\\\
\\epsilon_t &= \\sigma_t z_t, \\quad z_t \\sim N(0,1) \\\\
\\sigma_t^2 &= %s + %s\\, \\epsilon_{t-1}^2 + %s\\, \\sigma_{t-1}^2
\\end{aligned}
", mu, ar1, omega, alpha, beta)
}
```

在 Quarto 裡要讓它**渲染成數學式**（而非純文字），用 `results: asis` 的區塊把字串包在 `$$ … $$`：

````markdown
```{r}
#| results: asis
cat("$$", extract_garch_eq(fit), "$$")
```
````

注意：**係數名稱要對到實際模型**（先 `names(coef(fit))` 看有哪些，如 `ma1`、`omega`、`alpha1`、`beta1`、`shape`），
沒有的項就從模板拿掉；純 ARIMA 模型也可比照寫一個只填 `μ, φ_i, θ_j` 的版本。這比手寫更不易抄錯數字、也可重現。

## 不只金融：其他適合 SARIMA + GARCH 的領域

波動聚集（平時溫和、特定時期突然暴走且會持續一陣子）**不是金融獨有**。以下領域也常同時具有「均值的週期規律」
與「變異數的聚集」，適合 SARIMA + GARCH：**能源與公用事業**（電價、負載）、**氣象**（風速、降雨、溫度）、
**電商與網路流量**。

拿到一組**非金融的新資料**時，用兩個問題快速決策：

1. **資料有沒有「規律的週期，或前後互相影響的趨勢」？**
   - 有 → 用 **ARMA／SARIMA** 捕捉這個規律（**均值方程**）。
2. **資料有沒有「平時很溫和，但特定時期會突然暴走、且暴走會持續一陣子」的現象（即 ARCH 效應）？**
   - **有** → **務必加上 GARCH** 來建模／預測不確定性（波動度）。
   - **沒有**（波動範圍穩定）→ 不需要 GARCH，做到 **SARIMA 殘差為白噪音**即可收工。

判斷方法與金融分支相同：均值先用 SARIMA、殘差過一階 Ljung-Box；再對**殘差平方**用 `ACF` 與 `FinTS::ArchTest()`
檢定 ARCH——有才上 GARCH，沒有就停在白噪音殘差。

## 替代路線：Prophet（快速、強健的預測，偏預測導向 (B)）

**有 R 套件**：Meta 的 **`prophet`**（CRAN，`install.packages("prophet")`）；以及與 tidyverts 整合的
**`fable.prophet`**（`fable.prophet::prophet()`，可在 fable 工作流裡和 ARIMA／ETS 一起比較）。

**適用時機**：**假日／日曆效應、強季節性（可多重季節）、有缺漏或離群的商業數據**，
且**偏預測準確度／pipeline 穩健 (B) 取向**。Prophet 是可加性模型（趨勢 + 季節 + 假日），
對缺失與離群穩健、幾乎免手動定階、快速好上手，很適合零售銷量、網站流量等商業序列。

**不那麼適合**：需要**統計推論／可解釋性 (A)**、要交代參數顯著性與檢定時——Prophet 偏曲線配適與預測，
不提供 ARIMA 那種參數推論；這種情況仍以 ARIMA／SARIMA 為主。也因此，**波動度（ARCH/GARCH）不在 Prophet 範疇**，
若同時關心風險波動，Prophet 顧均值、變異數仍需另接 GARCH 思路。

基本用法（`prophet` 套件，資料需 `ds`（日期）＋ `y`（數值）兩欄）：

```r
library(prophet)
df <- data.frame(ds = dates, y = values)
m  <- prophet(df, yearly.seasonality = TRUE, weekly.seasonality = TRUE)
# 假日效應：用 holidays 參數，或 add_country_holidays(m, "TW") 加入國定假日
future   <- make_future_dataframe(m, periods = 30)
forecast <- predict(m, future)
prophet_plot_components(m, forecast)   # 分別看 趨勢 / 季節 / 假日 成分
```

屬 EDA 之後的建模銜接：EDA 階段先辨識「有沒有假日效應、強季節、缺漏」，再決定要不要走 Prophet。
