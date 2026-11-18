# RSDS 論文上の新規性と先行研究ギャップ（温熱環境 × 昼夜移動）

**2026年11月17日：公表済みの主要文献を照合。今回のPR #45の元GPS解析結果は原研究別に未完成であり、新規性の判断を肯定結果で上書きしない。**

## 主張の上限を先に修正

**「哺乳類は暑い日に夜行化する」は新発見ではない。** また「活動の総量と夜行性が同じ反応をしない」も既に実測されている。従って新しい時間ニッチの概念そのものを論文の主たるnoveltyにするのは弱い。

関連する根拠：

- [Gaynor et al.の人間活動変化に関する全球哺乳類比較, Nature Ecology & Evolution 2024](https://doi.org/10.1038/s41559-024-02363-2)。163種・102プロジェクトのカメラトラップで人間活動に伴う活動量と夜間比率の反応が地域／ギルドによって異なった。温度効果ではないが、人間活動・環境改変が時間ニッチと交絡し得ることを明示する必要がある。
- [Seeking temporal refugia to heat stress: increasing nocturnal activity despite predation risk, Proc. R. Soc. B 2024](https://doi.org/10.1098/rspb.2023.1587)。アイベックス47個体・2保護区の活動データで、高温後の夜間活動上昇と捕食リスクによる制約を報告。温度→夜間という機構の観測根拠は既にある。
- [A predominantly diurnal tropical mammal increases nocturnality in response to high temperatures, Biotropica 2023](https://doi.org/10.1111/btp.13271)。クチジロペッカリーのカメラトラップ＋微気候測定で気温との夜間検出比率の増加を報告。気温31.2°C以上で夜間活動率が気温1°Cあたり3.7%増えるとする閾値推定。**この活動検出比率とRSDSの速度中央値 `q` は別の単位／推論対象であり、効果量の大小を直接比較しない。**
- [Pallari et al., Temporal Responses to Warming: Do Wild Herbivores Trade Off Heat, Predators, and Humans?, Integrative Zoology 2026](https://doi.org/10.1111/1749-4877.70061)。地中海保護区のノロジカ・ダマジカのカメラトラップ（約12,400 camera-trap days、196地点）では両種の活動総量が暑い日に減る一方、夜行性の上昇はノロジカのみと報告。つまり**総量と時間配分の応答が一致しないこと、その種差の可能性は既に先行研究で実証されている**。
- [A guide to sampling design for GPS-based studies of animal societies, Methods in Ecology and Evolution](https://doi.org/10.1111/2041-210X.13999) と [Working with GNSS collar data. The importance of pre-analysis when setting the sampling interval, Ecological Informatics 2023](https://doi.org/10.1016/j.ecoinf.2023.102219)。GPS測位頻度と行動指標の関係は既知の方法学的課題。DAY/NIGHT各3区間の感度分析を「新しいGPS統計法」と売り込まない。

## 本研究の独自データ・生態学的検証で担えるギャップ

1. **同一の元公開GPSの速度成分比較**：カメラトラップ活動の夜間比率ではなく、元GPSをDAY/NIGHT同一ソーラー日に分けた`log1p(median observed step speed)`の夜間・日中**各速度変化と差 `q`**を同じ個体×年×月×気温・日長・先行水収支下で推定。
2. **原研究の独立性を支えるID・空間・日付の検証**：Utah／Wyoming／Subletteの3原研究DOIで夜間速度のTmax係数が正。元個体群の重複は未検証、WY内ミュールジカ・ビッグホーンの共通1°格子×年月は0。種間の結果を「同所条件での種応答差」と誤読しない。
3. **温熱関連がGPS昼夜の観測回数差だけで説明されるか**：DAY/NIGHT各3有効区間、同タグ／元日付／同ERA5に固定してもミュールジカの夜間気温係数の名目区間は3/3元研究で正。Utahの日中減速関連は、区間数をそろえると名目区間が0を含むため、測定への依存を正直に示す。
4. **今後の広域性**：少数研究から一般則を断定せず、Utah全4種・Wyoming全2種の原GPSを同一ネイティブモデルで再計算するPR #45。既存の旧Utah共通格子×年月のミュールジカ−ビッグホーン温度差は、水収支と日長を同時に入れると区間が0を含んだ。独立の共通条件の種間交互作用は未成立。

## アブストラクトの主張に使うべき表現

> Across selected original GPS telemetry studies of mule deer, modeled daily heat was consistently associated with an increase in observed nighttime displacement speed within tagged-animal calendar months. This directional association persisted when equal numbers of original GPS movement intervals were sampled during day and night; however, daytime speed responses were context-dependent, and broader interspecific differences were not established under shared geographic and seasonal opportunity.

これは **observational association** であり、暑熱回避・進化適応・夜行性の増加・直接活動時間・独立ランダム全球標本を意味しない。

**論文タイトル仮**：`Thermal associations of diel displacement differ across species and study contexts in original ungulate GPS telemetry`。元GPSの多研究による測定条件への頑健性と、その生態学的非一様性を中心にする。「universal nocturnalization under warming」は避ける。

## 未解決で、結論の新規性を左右する必須点

- **Utah全4種の同年月内原GPS結果は完了**：元全4種とWyoming2種の解析は成功し、元署名付き集約成果物の修正を[PR #46](https://github.com/zuizui0223/RSDS/pull/46)で再検証済み。NIGHT正5/6、q正3・負1・不確定2だが、独立原研究は2件。**Utahの同じ1°格子×同じ日付の温度差の種間検定**は[PR #47](https://github.com/zuizui0223/RSDS/pull/47)で実行中、結果が揃うまで種固有の熱反応の主張は保留。
- 独立したGPS標識個体群の空間／観測期間が重複しないことを公開原研究メタデータで監査する。
- **適切な比較**には、原GPS同じ1°セル×暦年×月で両種の観測機会が十分必要。Wyomingで両種の共有セル×年月0件は種間推論の強い制約。
- 1°C当たりの`q`効果量は`(1+vNight)/(1+vDay)`のモデル比率変化の尺度で概ね0.059–0.195%/°C（等3区間3研究）。カメラトラップによる活動の3.7%/°Cとは**まったく別の指標**であり、速度比の変化が生態学的に重要か、独立の解釈が必要。
- 太陽放射量→GPS移動方向の48外部原研究符号テストの結果 **NOT REPLICATED**、古いERA5 raw SHA12ソース **HOLD**は変更しない。

## 固定された解析と文書

- [3研究、原DAY/NIGHT各3区間・温熱関連の確定結果](../research/RSDS_THREE_MULE_EQUAL_GPS_DAY_NIGHT_FIRST_RESULT_20261116.json)
- [2元研究6ネイティブ種・旧元データ参照](../research/RSDS_UTAH_WYOMING_SIX_UNGULATE_NATIVE_HEAT_Q_REFERENCE_20261116.json)
- [PR #45・2元研究6種の新しい年月内原GPS](https://github.com/zuizui0223/RSDS/pull/45)
