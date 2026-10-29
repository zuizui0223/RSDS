# RSDS — 太陽・日射・時刻・季節の生態学

**RSDSを生態学研究の一次保管先とする。** 旧ODSPから、公開原データに基づく太陽方位・日射・日周移動・季節的再配分・等日長比較・マクロ気候解析を引き継ぐ。一般的なODSP状態空間予測アルゴリズムは本研究の中心にしない。

## 新しい鳥類の原GPS比較：ガン3種で日周移動ニッチと日長との関係が分離（2026-10-29）

**4公表Movebank元study（10登録分類群）の元GPS原ファイルを、行動結果を見る前にUUID・ファイルサイズ・MD5で固定し、4件すべてのチェックサム照合に成功。** そのうちガン3種の1元study（DOI 10.5441/001/1.fk899541）だけが、全登録分類群に共通する**120暦日・5年月ブロック**を満たした。カオジロガン62元タグの夜−昼GPS移動速度差qは**−0.18403**、マガン11タグ**−0.08391**、ハイイロガン10タグ**−0.06918**で、3種とも元個体bootstrap95%区間が負。一方、同じ元120日におけるq~研究共通参照日長の係数は**カオジロガン−0.02171/h、マガン−0.01581/h（不確定）、ハイイロガン+0.02066/h**。マガン−ハイイロガンの暦月内q種間差の日長係数は**−0.03594/h [−0.11046,−0.00392]**だが、5年月ブロックの事後・未補正探索的区間に限る。**時間ニッチの相対位置と季節的な方向が一つの生物指標に縮約できない**ことが、有蹄類・食肉類以外の新たな原鳥類studyでも観測された。

**3つの未成立原studyを削らない：** ガラパゴスのゾウガメは元ファイルMD5一致でも登録学名の厳密行照合0件のため**学名schema HOLD**（動物不在ではない）。カモではマガモ21タグの元q負が推定できたが2種共通適格日0、アフリカのハゲワシは昼夜GPS区間支持が不十分。3種に共通する元暦日は同所共存・実際の日射経験の同一性を意味せず、qも真の活動時間ではない。既存48原studyの日射外部負方向非再現と旧ERA5 raw JSON12コンテキストHOLDは変更しない。[結果の原4study×10 native taxa JSON](research/RSDS_FOUR_NEW_MULTINATIVE_GPS_SHARED_CALENDAR_V1_FIRST_RESULT_20261029.json)／[SHA固定したガン類120日・全HOLD検証](research/RSDS_FOUR_MULTINATIVE_ORIGINAL_GPS_ECOLOGICAL_VALIDATION_V1_20261029.json)／[生態学的解釈](docs/RSDS_GOOSE_THREE_NATIVE_SHARED_DATE_TEMPORAL_NICHE_20261029.md)／[原GPS・検証CI](https://github.com/zuizui0223/RSDS/actions/runs/40860577748)。


## 新しい生態学的検証：日長と時間ニッチの関係を同じ暦月内で比較（2026-10-29）

**Illinois 3種の169共通暦日、Wyoming 2種の193共通暦日**で、原Movebank GPSのMD5・元q／日長係数・種間共通日付を再現後、43°N共通天文日長と種間q差を同じ**暦年×暦月内で中心化**した。6つの固定種間×GPS時計比較を全て報告し、年月ブロック2000回の記述的bootstrap95%区間が**2比較で従来と同符号、4比較でゼロを含んだ**。支持されたのはIllinoisボブキャット−オジロジカ **−0.08466 [−0.17655,−0.00108]**、Wyomingミュールジカ−ビッグホーン（元GPS時刻） **+0.01027 [+0.00008,+0.02290]**。他4比較の**点推定はゼロではない**が、月ブロックの不確実性で符号が未確定。年月ブロックは各研究8と少数で、日長と日付・季節資源・気象はまだ分離できない。**これは2原研究内の事後的な観測的結果で、光周期の因果効果・世界的普遍性ではない。** [元結果JSON](research/RSDS_SHARED_MONTH_ADJUSTED_PHOTOPERIOD_TWO_STUDY_V1_FIRST_RESULT_20261029.json) ／ [SHA固定の結論検証](research/RSDS_SHARED_MONTH_PHOTOPERIOD_VALIDATED_ECOLOGY_20261029.json) ／ [生態学的解釈](docs/RSDS_ILLINOIS_WYOMING_WITHIN_MONTH_PHOTOPERIOD_ECOLOGY_20261029.md)。旧48原研究の日射外部再現不成立とERA5旧JSON12元HOLDは不変。

## 新しい原GPS比較：同じ観測日でも3種の時間利用が異なる（2026-10-29）

**Illinoisの公開Movebank元GPS研究（DOI 10.5441/001/1.649）で、コヨーテ・ボブキャット・オジロジカが同時に観測された169暦日／9年月ブロックを抽出した。** 元ファイルUUID・bytes・MD5、3種すべての既報qと日長係数を**新しい同暦日比較の前に数値再現**。適格元タグはコヨーテ8、ボブキャット5、オジロジカ12。

夜−昼の観測GPS直線変位速度差qの平均は**コヨーテ +0.26852 ＞ ボブキャット +0.09339 ＞ オジロジカ +0.00245**。3つの種間差すべてで元タグbootstrap95%区間は正。一方、共通43°N天文日長に対するqの係数は**コヨーテ +0.03744（正）・ボブキャット −0.02343（負）・オジロジカ −0.00037（不確定）**。**同じ元研究・共通暦日でも、日周移動の相対位置と季節的関連は別々の特徴として変わる**という観測的結果が得られた。

[元3種・169共通暦日の結果](research/RSDS_ILLINOIS_THREE_NATIVE_SHARED_DATE_V1_FIRST_RESULT_20261029.json) ／ [生態学的解釈](docs/RSDS_ILLINOIS_COYOTE_BOBCAT_DEER_SHARED_TEMPORAL_ECOLOGY_20261029.md) ／ [原GPS取得・既報値再現のCI](https://github.com/zuizui0223/RSDS/actions/runs/40841033373)。**制約：** 独立な原フィールド研究は**1件**、ボブキャットは5タグ、元GPS移動差qは真の活動時間割合でない。暦日が同じでも微小生息環境・採餌・捕食行動・日中の観測区間は一致しない。事後的な選択と多重比較未補正の名目区間であり、直接的な捕食・競争や日長因果効果ではない。**別のWashington原GPS研究（DOI 10.5441/001/1.gm93267b）は両種のq適格タグが0で地域再現はHOLD**。48元研究の日射方向の外部負方向非再現、ERA5旧12コンテキストのHOLDはいずれも不変。

## 最新の独立原GPS比較：日周移動の順位は2地域で共通、日長応答の種間差は共通しない（2026-10-29）

**Utahの別公表Movebank原GPS（DOI 10.5441/001/1.711）で、Wyomingと同じ2種（ミュールジカ・ビッグホーン）を事前固定し、元研究内の199共通暦日・2年・10／42元タグを検証。** 原GPS提供元MD5、Utah全4分類群の以前のq/日長結果、選択した2種の3測位時刻系列の係数を新比較前に再現した。Wyomingは193共通暦日・108／41元タグで、両原研究とも同じ43°N天文日長を採用。

**新しい肯定的結論：** 共通暦日の観測NIGHT−DAY GPS直線変位速度差qの*ミュールジカ−ビッグホーン*の種間順位は、元系列・UTC偶数時2h・奇数時2hの**2研究×3系列すべてで正の元タグbootstrap95%区間**（Wyoming元系列 **+0.01702 [0.01416,0.01987]**、Utah元系列 **+0.03767 [0.02158,0.05716]**）。一方、q〜日長傾きの**ミュール正・ビッグホーン負**というWyomingの対照（種間差 **+0.01235 [0.00982,0.01496]**）はUtahでは支持されず、Utahの2種の元系列傾きはともに正で、種間差 **−0.02720 [−0.07008,+0.01320]** はゼロを含む。**日周ニッチの相対的な位置と季節的な日長関連は、異なる範囲で一般化する可能性**がある。

[Utah・Wyomingの2原study×共通暦日の再現検証](research/RSDS_UTAH_WYOMING_TWO_STUDY_SHARED_DATE_PHOTOPERIOD_V1_RESULT_20261029.json) ／ [生態学的解釈](docs/RSDS_UTAH_WYOMING_DIEL_RANK_VS_PHOTOPERIOD_SLOPE_20261029.md) ／ [Utah元GPS first result](research/RSDS_UTAH_COOCCURRENT_DATES_PHOTOPERIOD_V1_FIRST_RESULT_20261029.json)。**制約：** Utahのミュールジカは10元タグのみ、元GPSのUTC測位時刻に強い偏り。種の選定はWyomingの先行結果を見た後であり、完全なブラインド外部再現ではない。qは真の活動時間ではなくGPS変位速度代理指標で、日長の因果作用や世界一般則は未同定。旧ERA5 12HOLD・48元研究の日射外部符号非再現はそのまま。

## 最新の肯定的生態学的結果：同じ調査・同じ暦日でも2種の日長関係が反対（2026-10-29）

**Wyomingの1つの公開Movebank原GPS研究（DOI 10.5441/001/1.703）から、ミュールジカとビッグホーンがともに観測された193暦日・2年・8年月ブロックを元データから復元した。** 元GPSのMD5、両種の過去の全データ／3系列共通日付のq・日長係数を**新比較前に再現**。同一の**43.0°N基準天文日長**を使い、元の1時間UTC bin系列・偶数UTC時2時間系列・奇数UTC時2時間系列のすべてで、**ミュールジカのNIGHT−DAY観測GPS速度差q~日長係数は正、ビッグホーンは負**だった。共通暦日で適格な原タグは**108頭対41頭**。元系列の種間日長係数差は**+0.01235／日長1時間（元タグbootstrap95% [+0.00982,+0.01496]）**。さらに193暦日の種別平均qを同日ごとに引いた比較でも、**+0.01167（8年月ブロックbootstrap [+0.00763,+0.02029]）**の正の傾きとなり、偶数・奇数2時間系列でも正。**種別の観測年・季節が全く異なるだけ**では説明しにくい、研究内の時間ニッチの観測的異質性が残った。

[元Wyomingの検証済み結果](research/RSDS_WYOMING_COOCCURRENT_DATES_PHOTOPERIOD_V1_FIRST_RESULT_20261029.json) ／ [SHA固定の定量監査](research/RSDS_WYOMING_SHARED_SEASON_DIVERGENT_PHOTOPERIOD_VALIDATION_V1_20261029.json) ／ [科学的な解釈と限界](docs/RSDS_WYOMING_SHARED_CALENDAR_PHOTOPERIOD_ECOLOGY_20261029.md) ／ [統合論文草稿](manuscript/RSDS_EIGHT_STUDY_TIME_SOLAR_MOVEMENT_MACROECOLOGY_v1_20261026.md)。

**推論の上限：** 1つの原研究の事後感度検証であり、2種の生息場所・元タグ・測位機器と時別気象条件が完全一致したわけではない。両種に共通する193暦日のうち各元タグが観測した日は部分集合。日長は事前固定の43°Nで標準化した天文時計であり、動物が受けた日射量ではない。qは実際の活動時間割合ではなくGPS直線変位速度代理指標。適応・因果的日射反応・世界一般性を実証したことにはならない。旧48元studyの日射方向外部非再現と旧ERA5 JSON12保留も不変。

## 最新の生態学的結果：観測時刻をそろえても残る時間ニッチの違い（2026-10-28）

**元Movebank 8公表研究・12調査×分類群の原GPSを、元の全観測日と同一元個体・年・日付に厳密照合した3つのUTC測位時刻系列で比較した。** [実測生態学的係数](research/RSDS_TRIPLE_DATE_ECOLOGICAL_Q_DAYLENGTH_INVARIANCE_V1_FIRST_RESULT_20261028.json)／[元全データとの厳密な比較](research/RSDS_TRIPLE_DATE_VS_ORIGINAL_ECOLOGICAL_TIME_SIGN_AUDIT_20261028.json)／[解釈](docs/RSDS_SAME_TAG_DATE_TEMPORAL_NICHE_SIGN_ROBUSTNESS_20261028.md)。**元の全データと同一個体・同じ日付の3系列のすべてでNIGHT−DAYの観測GPS移動速度差qの符号が一致したのは5/8比較可能群**（夜側の相対GPS移動大：ボブキャット・英国アカギツネ・コアラ。昼側：UtahとWyomingのビッグホーン）。**q~天文日長傾きの符号まで同じ条件で正として支持されたのは2/8群**（Wyomingミュールジカ、カリフォルニア・ボブキャット）。残りの結果には共通元タグの日数不足と名目bootstrap区間の不確定性がある。特にUtahミュールジカでは全適格タグ1,723のうち3系列共通日付タグは11のみで、結果不確定を生物学的否定に置き換えない。これは**GPS直線変位の時間ニッチ代理指標の符号に限る記述的な頑健性**であり、実活動時間、日射・日長の因果作用、世界一般則は未同定。旧48原日射研究外部非再現と旧ERA5バージョンHOLDは維持する。

## ERA5気象の実際の数値版の再現性（2026-10-28）

**新しい独立検証：** 動物GPSではない**固定2地点のERA5全時別気象**を各2回取得すると、**2/2地点で生JSON SHA256が変化した一方、UTCの全8760/8784時間・短波放射・気温・風速・風向の数値版SHA256は同一**、返却モデル格子も一致した。つまりファイルのバイト列のSHA不一致だけでは「実際の気象数値が変更された」と言えない。[検証済みベンチマーク結果](research/RSDS_ERA5_CANONICAL_BENCHMARK_FIRST_RESULT_20261028.json)／[仕様と科学的境界](docs/RSDS_ERA5_PHYSICAL_VS_JSON_VERSION_AUDIT_20261028.md)。従来の**元GPS×ERA5 12コンテキストの旧JSON SHA不一致はHOLDのまま**で、当時の数値配列が未保存なので旧と現在の物理数値が等しかったかは判定不能。以前の外部負方向日射効果の非再現も変更しない。

**さらに実際の原GPS解析の12環境コンテキストまで検証（2026-10-28）：** Utah 4分類群・Wyoming 2分類群の元Movebank GPSファイルMD5、q・日長係数・適格個体数を**ERA5取得前に再現**。以前に固定した**12コンテキスト（重複除去で10グリッド年）**をそれぞれ2回照会した結果、**生JSON SHAは12/12で不一致、時別の気象4物理変数・UTC時刻・単位・モデル格子は12/12で一致（物理数値が異なる時点0）**。旧12コンテキストの旧生JSON SHAとの一致は0だが、旧物理時系列が保存されていないため**過去の物理値の再現性は不明**で、旧ERA5_VERSION_HOLDは変更しない。[12元GPSコンテキストの計測結果JSON](research/RSDS_ERA5_ORIGINAL_SOURCES_CANONICAL_FIRST_RESULT_20261028.json)／[解釈](docs/RSDS_ORIGINAL_GPS_ERA5_PHYSICAL_CONTENT_12_CONTEXTS_20261028.md)。[GitHub Actions検証](https://github.com/zuizui0223/RSDS/actions/runs/40534683501)。

## 日射残差と移動方向：新ERA5物理版で元GPS統計係数を再現（2026-10-28）

**元GPSの科学的係数が再計算可能になった：** Utah/Wyomingの全6分類群、元の12環境コンテキスト×3軸（昼夜GPS速度差q、太陽相対方位cosとsinへの短波放射残差の関連）を、[新たに正規化・固定したERA5物理数値版](research/RSDS_ERA5_ORIGINAL_SOURCES_CANONICAL_FIRST_RESULT_20261028.json)で再解析。**以前保存した36/36の条件付き係数と最大絶対差2.1×10⁻¹⁷以内、6群×3軸の18/18群平均も最大差1.4×10⁻¹⁷以内で再現。** [数値比較JSON](research/RSDS_CANONICAL_ERA5_V2_OLD_V0_EFFECT_REPLAY_AND_TAG_PAIRS_20261028.json)／[新結果JSON](research/RSDS_CANONICAL_ERA5_SAME_TAG_ECOLOGY_V2_FIRST_RESULT_20261028.json)／[解釈](docs/RSDS_CANONICAL_ERA5_SAME_TAG_ECOLOGY_V2_20261028.md)。同じ元タグの2コンテキストで4/6群が比較可能となり、UtahワピチとWyomingミュールジカのq係数差、Wyomingビッグホーンの太陽直交sin係数差に未補正の正の区間がある。ただし元GPS観測支持順位2−順位1であって**高日射・低日射の因果差ではない**。**旧ERA5の時別物理数値が保存されていないこと・旧12生JSON SHAのHOLD・48元研究の負方向外部非再現判定はそのまま**。同じ2公表研究の事後再解析であり、新しい独立原調査での生物学的再現と数えない。

## GPS観測時刻の偏り：8原調査からの新結果（2026-10-28）

**元原GPSを検証した8研究×12分類群の1時間bin・UTC偶数/奇数2時間系列の比較が成功。** [全study結果](research/RSDS_EIGHT_STUDY_GPS_CADENCE_MATCHED_DAY_FIRST_RESULT_20261028.json)／[UTC時刻偏りの固定診断](research/RSDS_EIGHT_STUDY_GPS_CADENCE_PHASE_ALIAS_DIAGNOSTIC_20261028.json)／[生態学上の解釈](docs/RSDS_GPS_SAMPLING_CADENCE_SAME_TAG_ECOLOGY_20261028.md)。Utahの4分類群ではGPS点の**87.1〜95.6%が奇数UTC時**に偏り、元の「1時間binから最大1件」の系列が実際の1時間ごとの受信を意味しないことが判明した。Wyomingのミュールジカ・ビッグホーンはUTC偶奇がほぼ均等だが、偶数時・奇数時の元系列との差でq（観測NIGHT−DAYのGPS速度差）の**記述的符号が双方の種で反転**した。両位相の採用個体・日付が同じとは限らず、同一個体の行動逆転ではない。野生シチメンチョウは2時間側でq適格日が0、ダイサギは元ファイル容量上限HOLD。現在は3系列すべてで同じ元タグ・日付に限定する追加確認へ進める。**旧ERA5版HOLDと48元研究の日射負方向外部非再現は不変。**

**厳密な同じ標識個体・年・日付の検証も完了（2026-10-28）：** [追加原GPS 3系列検証の第一次結果](research/RSDS_GPS_THREE_PHASE_EXACT_TAG_DATE_V2_FIRST_RESULT_20261028.json)／[固定SHA実証整理](research/RSDS_GPS_TRIPLE_DATE_CADENCE_ECOLOGICAL_INTERPRETATION_V2_20261028.json)／[科学的解釈](docs/RSDS_GPS_EXACT_THREE_PHASE_TAG_DATES_ECOLOGY_20261028.md)。8原study×12 study/taxon行を保持し、元系列・偶数UTC時・奇数UTC時の**3系列全てに共通する元タグ・元暦年・元天文日**を使って確認した結果、**8群で元タグ5個体以上の完全対応が成立**。Wyomingミュールジカ**113元タグでq偶数−奇数 +0.005393 [0.003574,0.007132]**、ビッグホーン**41元タグで+0.013339 [0.011617,0.014897]**と、記述的元個体bootstrap95%区間がゼロを含まなかった。**異なる個体／日付を採用したことだけでは説明できないGPS測位時刻・移動区間の系統的差**である。ただし偶数時と奇数時で実際に観測したGPS移動区間は異なるため、動物の行動や太陽コンパスがUTC偶奇で変化した証明ではない。元study単位は2件（Wyoming+Utahの陽性3群）に過ぎない。古いERA5 HOLDや日射効果の外部負方向非再現は不変。

## 現在の中心質問：時間利用と方向利用の組織化

> **太陽・季節・環境の変化に対して、生物の時間利用と移動方向はどのように組織化され、その違いは種や地域を越えて一般化するのか？**

日長・太陽方位・太陽高度・放射・気温・風を「環境状態」とし、観測された**昼夜の相対移動速度（時間ニッチのGPS代理量）**と、**太陽相対GPS移動方向（移動ニッチ）**を分ける。両者が同じ原調査で比較可能か、さらに別の原調査にも残るかを中心に検証する。撮影時刻は検出の観測過程を明示し、真の活動時間と区別する。

**2026-10-26 最新のマクロ生態学的拡張（原GPS 8調査）：** 以前のUtah/Wyomingの2元調査と、新たにGPSイベントを開く前に原ソースを固定した**6元Movebank調査**を統合。**8元study DOI・12 study×taxon行・10学名**の比較表と[再現可能なJSON](research/RSDS_8_STUDY_TIME_SOLAR_HEADING_MACRO_ATLAS_V1_20261026.json)を確定。[研究内容と非支持結果の詳細](docs/RSDS_EIGHT_STUDY_TEMPORAL_DIRECTION_MACROECOLOGY_20261026.md)／[英語論文草稿](manuscript/RSDS_EIGHT_STUDY_TIME_SOLAR_MOVEMENT_MACROECOLOGY_v1_20261026.md)。NIGHT−DAYの観測GPS最短移動速度差qは**11/12行で推定可能（負7・正4）**、日長とqの個体年内係数は**5正・1負・5未確定（1元ファイルHOLD）**。生の太陽相対cosは**10/12行推定可能（2正・2負・6が0を含む）**。同一元タグの2観測コンテキストは**9/12行・延べ288タグペア**。世界代表の12独立生態系ではなく、元タグ間多重比較・機器観測条件も未調整。**新6調査のうちダイサギの原GPSファイルは事前200MB上限を超えHOLD**、スナドリネコの太陽方向の支持個体は3で未推定。古い48原調査の日射外部再現失敗およびERA5 JSON差の12保留は変更していない。

**2026-10-26 元研究1件＝1単位の再集計（8研究の階層感度）：** [source-equal生態学ノート](docs/RSDS_EIGHT_STUDY_SOURCE_EQUAL_MACRO_SENSITIVITY_20261026.md)を追加。Utahの4分類群、Wyomingの2分類群をそれぞれ1原調査にまとめると、q=観測されたNIGHT−DAYのGPSログ変位速度差の平均は、**調査×分類群11行の均等平均 +0.00480**から、**元study7研究の均等平均 +0.03464**に変化した。研究均等の記述的なbootstrap95%区間 **[−0.03414,+0.10346]** はゼロを含み、全球で夜行性が優勢という結論ではない。旧6群のjoint solar covarianceは原タグ区間を欠くため新しいCIを捏造せず**点推定に限定**。元資料・外部再現失敗・HOLDは不変。[固定計算結果](research/RSDS_EIGHT_STUDY_SOURCE_EQUAL_MACRO_SENSITIVITY_FIRST_RESULT_20261026.json)。

**2026-10-26 新しい原GPS統合実証：** [RSDSの時間×方向マクロ生態学の結果](docs/RSDS_JOINT_SOLAR_TEMPORAL_MOVEMENT_MACROECOLOGY_20261026.md)を公開。[実測の固定JSON](research/RSDS_CROSS_AXIS_GENERALITY_FIRST_RESULT_20261026.json)、[再現コード](scripts/analyze_rsds_cross_axis_generality_v0.py)、[論文草稿](manuscript/RSDS_SOLAR_TIME_DIRECTION_MACROECOLOGY_v0_20261026.md)および[独立CI](https://github.com/zuizui0223/RSDS/actions/runs/40262509280)を保持。Utah/Wyomingの原GPSでミュールジカのq−日長傾きは**+0.005619/+0.005538**と両方正である一方、太陽相対の生方向平均は**+0.023393/−0.048758**と逆符号。6原study×taxon行、9学名の2ソース方位再現監査、48原ソースの日射効果の対応をDOI・学名で機械的に照合した。

**2026-10-26 同一個体による時間・方向・日射の新しい実証：** [同じ標識個体・観測日の原GPSと時別ERA5の研究](docs/RSDS_SAME_ANIMAL_SOLAR_IRRADIANCE_MOVEMENT_PROTOCOL_20261026.md)を実行。Utah 4分類群、Wyoming 2分類群の原GPSのプロバイダMD5と、全6群の既存の昼夜移動・日長係数を完全再現したうえで、**元タグ個体434の時間・方向の両関連**を計測した。[原結果JSON](research/RSDS_SAME_ANIMAL_SOLAR_WEATHER_V0_FIRST_RESULTS_20261026.json) ／ [科学的比較結果](research/RSDS_SAME_ANIMAL_SOLAR_WEATHER_V0_ECOLOGICAL_SYNTHESIS_20261026.json) ／ [英語原稿](manuscript/RSDS_SAME_ANIMAL_IRRADIANCE_TIME_DIRECTION_v0_20261026.md)。日射の天文幾何・1時間太陽時刻・風・気温を調整すると、全6群で**昼夜相対GPS変位と日射残差との係数の個体bootstrap区間が0を含む**一方、Utahビッグホーンの太陽相対方向cos係数は **+0.01686 [0.00636,0.02647]**、sin係数は **+0.02917 [0.01765,0.04126]** と元調査内で正の記述的関連を示した。Wyomingビッグホーンでは元2グリッド年のsin成分係数が逆符号となった。

**2026-10-26 ERA5に依存しない、同一タグ個体の新しい比較（原GPS v2）：** [検証済みの生態学的結果](docs/RSDS_SAME_TAG_GPS_GEOMETRY_ECOLOGY_20261026.md)を追加。UtahとWyomingの2原調査・6調査×分類群で、**元の同じタグ個体の観測コンテキスト2つを対応づけた延べ269個体ペア**を確保した。観測されたNIGHT−DAY GPS速度差 q は**5/6群で個体bootstrap95%区間が0を含まない**（正2、負3）。生の太陽相対移動cos差は**2/6群**で区間が0を含まないが、その2群では地理的移動方向marginalも異なり、**太陽方向とのjoint circular covariance差は6/6群で区間が0を含んだ**。未補正の記述的区間で、太陽コンパスや同じ環境勾配の因果的応答は未証明。元GPSだけの解析であり、以前のERA5再取得JSON不一致**12/12の保留**もそのまま保持。[原結果JSON](research/RSDS_SAME_TAG_GPS_GEOMETRY_V2_FIRST_RESULT_20261026.json) ／ [統合結果](research/RSDS_SAME_TAG_GPS_GEOMETRY_V2_ECOLOGY_SYNTHESIS_20261026.json) ／ [論文原稿](manuscript/RSDS_SAME_TAG_GPS_TEMPORAL_SOLAR_GEOMETRY_v2_20261026.md)。

**主張の境界：** 従来の48原日射調査とUtah/Wyomingの2調査では元DOIが重ならなかったが、今回の再解析により**後者2調査内**では同一個体・観測日への環境結合を実現した。これで世界一般の因果的「太陽方向ニッチ」が識別されたわけではなく、源元個体・地域・原GPSの移動経路・選択したグリッド年と気象履歴には交絡が残る。事後探索の複数群比較であり、旧48原調査の事前指定外部再現判定 **NOT REPLICATED** は一切変更しない。

**2026-10-26移管：** [移転元の固定コミット](https://github.com/zuizui0223/odsp/commit/d58fb759416fc458f69ef61ea79b7b1d4164f2c5) の太陽・季節生態学に関する結果、契約、スクリプト、テスト、草稿、台帳、歴史的CI仕様を保存する。[固定390ファイルの移管マニフェスト](provenance/ODSP_SOLAR_SEASON_MIGRATION_MANIFEST_20261026.json)と[移管照合記録](provenance/RSDS_TRANSFER_VERIFICATION_20261026.json)を参照。照合が全件成立するまで「完了」と見なさない。

## 中心となる生態学的問い

1. **同じ日長でも、動物の昼夜移動は同じか。** 上昇・下降日長枝を個体内で比較し、気象・植生・季節的環境の非同期性を検証する。
2. **動物の移動方向は、太陽方位に対してどう分布するか。** 種・個体・地域・季節の違いを原GPSで比較する。
3. **日射量・温度・風向・太陽位置の関係は、動物の方向選択とどう結びつくか。** 太陽高度・細かい太陽時刻・移動経路を分離する新しい実データ比較を進める。

## 現在までのソース検証済み重要結果

| テーマ | 原データからの肯定的結果・異質性 | 主張の境界 |
|---|---|---|
| 同日長の昼夜移動 | Movebankの10原ソースコホート・254適格個体で昼夜相対GPS変位の季節枝差は一様でない。ワピチ +0.045、アメリカクロクマ −0.077など。 | サンプルは非ランダムで個体の真の全活動時間配分とは同一でない。 |
| 独立2原調査での比較 | ミュールジカの日長−昼夜相対変位の個体内・年内係数はUtah +0.00562、Wyoming +0.00554。 | 2元調査での記述的再現。日長の因果作用と一般化は未確定。 |
| 等日長の気候状態 | 固定した21陸上地点×3年で有効62/63のERA5ペアのうち29ペアで気温と表層土壌水分の変化方向が反対。 | 任意に固定した環境再解析点。全球代表、動物行動の根拠ではない。 |
| 太陽相対GPS方向 | 原GPS40ソース・41分類群コホートで太陽との角度関係に異質性。 | 太陽コンパス利用・日射回避の直接証明ではない。 |
| 日射量と移動方向 | 原GPS48ソース中20で日射係数推定可。初期10ソース平均 −0.01966／+100 W/m²、外部10平均 −0.00336。 | **事前指定した負方向の外部再現基準は不成立。**28非推定をゼロにしない。 |
| 同一元GPS個体の時間×方向×日射（探索） | 元2調査・全6分類群、計434個体で両軸を同じ観測日・ERA5時刻へ結合。Utahビッグホーンのcos・sin方向係数は正の個体bootstrap区間。 | 日射と昼夜相対移動係数は6群すべてで区間が0を含み、普遍的な光忌避・環境媒介・競争効果は未確定。 |
| マクロ気象勾配（探索） | 20有効原ソースの係数と気温rho +0.448、風流向−太陽一致度 +0.499、緯度近似 −0.564、平均日射 +0.075。 | 事後探索・非ランダムソース、個体内因果係数ではない。 |
| 原カメラ写真時刻 | Serengeti等でほぼ同日長の異季節枝でも記録写真時刻分布の分類群差。 | 検出感度と稼働時間の補正不足のため、真の行動時計とは区別。 |

**完全な結論と失敗も含む境界：** [RSDS総合研究結論](docs/RSDS_ECOLOGICAL_SYNTHESIS_20261026.md)

## 内容と再現

- [docs/](docs/)：原解析レポート。 [manuscript/](manuscript/)：研究草稿。
- [research/](research/)：元の研究契約、肯定・不成立・HOLDの完全結果、Providerチェックサム、[Movebank全ソース台帳](research/seasonal_database_catalog/2026-10-19/movebank/)。
- [scripts/](scripts/)：原データに再接続する研究用スクリプト。 [tests/](tests/)：既存の固定事実テスト。
- [archive/odsp_workflows/](archive/odsp_workflows/)：旧ODSPで利用された歴史的CI定義。RSDSでは直接起動しない。
- [provenance/](provenance/)：移管元commit、各ファイルのblob SHA、RSDS側の一致判定。

ローカルのマクロ生態学再計算は Python 3.11 以降で、次のコマンドを使用する。

    python -m unittest discover -s tests -p 'test_solar_irradiance_direction_macro_gradient_v0.py' -v
    python -m scripts.analyze_solar_irradiance_direction_macro_gradient_v0 --out outputs/solar_irradiance_direction_macro_gradient.json

**再利用と研究倫理：** 元の動物GPS座標・ID・正確なイベント時刻・原写真をRSDSに再配布しない。公開元DOI、プロバイダMD5、元ソースの利用許諾は保持。公開Movebankカタログの取得ページ圧縮ファイル（raw_pages.tar.gz）は二次的アーカイブであり、バイナリではなくそのCSV/JSON台帳とmanifestを移管する。ODSPの旧結果は改ざんせず履歴として残すが、**日射・太陽・時間・季節の生態学的解釈と新規発展の正本はRSDS**とする。
