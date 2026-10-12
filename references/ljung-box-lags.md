# Ljung–Box 檢定的滯後期數 m：Tsay 的實際用法與本 skill 的預設

> 本檔記錄「m 要取多少」的**查證結果**。查證來源：Tsay《Analysis of Financial Time Series》講義
> Ch.2（lec3，linear time series）與 Ch.3（lec4，conditional heteroscedastic models），以及 `FinTS::AutocorTest` 文件對 Tsay 的引用。
> 流行說法「Tsay 用 m = 10 做主要檢定；GARCH 後對殘差平方用 m = 20」**部分正確、部分不準確**，詳見下表。

## 1. 講義裡 Tsay 實際用過的 m

| 情境 | 講義 | m |
|---|---|---|
| ARMA 殘差（J&J 季資料 airline 模型） | lec3 第 41 頁 | **12**（Q(12) = 10.0，p = 0.44） |
| 報酬序列本身（月報酬，檢查有無自相關） | lec4（Intel 例） | **12**（Q(12) = 18.26，p = 0.11） |
| ARCH 效應檢定（原始報酬平方 `a_t²`） | lec4（Intel 例） | **12**（`ArchTest(x, lags = 12)`；Q(12) = 89.85） |
| ARCH(3) 後，標準化殘差 `ã_t` 與其平方 `ã_t²` | lec4 | **10**（Q(10) = 12.64；Q(10) = 14.75） |
| GARCH(1,1) 後，`ã_t` 與 `ã_t²` | lec4 | **10 與 20 都報**（`ã_t`：Q(10) = 11.22、Q(20) = 24.29；`ã_t²`：Q(10) = 9.92、Q(20) = 16.74） |
| EGARCH 後，`ã_t` 與 `ã_t²` | lec4 | **10 與 20 都報**（`ã_t`：6.65／20.52；`ã_t²`：3.37／11.15） |

## 2. 對流行說法的判讀

| 說法 | 結論 | 依據 |
|---|---|---|
| 「Tsay 用 m = 10 作為主要檢定」 | **部分正確** | 在波動度模型的診斷（ARCH／GARCH／EGARCH 後）確實以 Q(10) 為首；但**月、季資料的均值檢定他用 m = 12**（對應一個年度週期的倍數） |
| 「GARCH 後對**殘差平方**用 m = 20」 | **不準確** | Tsay 的做法是 `ã_t` **與** `ã_t²` **各自**都報 Q(10) 和 Q(20)，並非「平方用 20、其他用 10」。20 是**補充較長期的相依**，不是只給平方用 |
| 「對**殘差平方**檢定」 | **要修正對象** | 該對**標準化殘差** `ã_t = a_t / σ̂_t` 取平方，不是原始殘差的平方。Tsay：`ã_t` 的 Q 檢**均值方程**是否足夠，`ã_t²` 的 Q 檢**波動度方程**是否足夠 |
| 理論建議 | 補充 | Tsay 書上指出模擬研究顯示 **m ≈ ln(T)** 的檢定力較佳；有季節時要改為包含季節週期的倍數（見 `FinTS::AutocorTest` 文件對 Tsay 的引用）|

> `m ≈ ln(T)` 只是一個起點：T = 2467 時 ≈ 7.8，T = 250 時 ≈ 5.5。Tsay 自己的例子多取 10 或 12，所以 10 是**慣例**，不是定理。

## 3. 本 skill 的預設（依上述查證）

| 檢定對象 | 預設 m |
|---|---|
| **一階矩**：ARMA／ARIMAX／迴歸殘差、報酬序列的 Ljung–Box | **m = 10**；有季節性時，m 至少涵蓋一個完整季節週期（月資料 12、季資料 12 為 Tsay 慣例）|
| **GARCH 配適後**：標準化殘差 `ã_t` **及其平方** `ã_t²` | **各報 Q(10) 與 Q(20)**（共 4 個統計量），先看 10，再用 20 確認較長期 |
| **ARCH 效應檢定**（GARCH 之前，檢原始殘差／報酬的平方） | **m = 12**（Tsay 的 `ArchTest(lags = 12)` 與 Q(12)）；資料很短時可降低 |

任何情況都必須 **m > 已估計的參數個數**，且報告要寫明所用的 m（假說寫法見 SKILL.md）。

```r
# GARCH 配適後：標準化殘差 ã_t 與 ã_t²，各檢 m = 10、20
z <- as.numeric(residuals(fit, standardize = TRUE))     # rugarch：標準化殘差
for (m in c(10, 20)) {
  q1 <- Box.test(z,   lag = m, type = "Ljung-Box")      # 均值方程是否足夠
  q2 <- Box.test(z^2, lag = m, type = "Ljung-Box")      # 波動度方程是否足夠
  cat(sprintf("m = %2d | ã_t: Q = %.2f (p = %.3f) | ã_t^2: Q = %.2f (p = %.3f)\n",
              m, q1$statistic, q1$p.value, q2$statistic, q2$p.value))
}
```

## 4. 尚未查清的細節

- **自由度**：用講義 p 值反推所用的 df：
  - GARCH(1,1) 例：`ã_t` 的 df = m（10、20），`ã_t²` 的 df ≈ m − 2（8、18）。
  - ARCH(3) 例：兩者 df 都是 10，**沒有調整**。
  - EGARCH 例：`ã_t` 的 df ≈ m − 1（9、19），`ã_t²` 的 df ≈ m − 2（8、18）。
  - **講義文字未說明規則，以上是從 p 值反推，不是 Tsay 的明文**；Tsay 並非每例都調整。上方範例程式未設 `fitdf`（df = m）。
- `rugarch` 的 **Weighted Ljung–Box**（Fisher & Gallagher 2012）要用 **`show(fit)`**（或直接印出 `fit`）才會輸出；`summary(fit)` 只回傳一般摘要表，**不含**這些診斷。
  它的落後期數由 rugarch 依 p+q 自訂（GARCH(1,1) 實測：標準化殘差 Lag 1、2、5；平方 Lag 1、5、9），**不是 Q(10)、Q(20)**，與上面手動檢定**不必一致**；兩種都可以報，但要註明哪一種。

## 參考

- Tsay, R.S. *Analysis of Financial Time Series*, 3rd ed.；課程講義 lec3（Ch.2.7–2.11）、lec4（Ch.3）。
- `FinTS::AutocorTest` 文件（引述 Tsay：m ≈ ln(T)；季節資料需修正；df = m − g）。
- Fisher, T.J. & Gallagher, C.M. (2012). New weighted portmanteau statistics for time series goodness of fit testing. *JASA*.
