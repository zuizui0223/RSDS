# RSDS — 太陽・日射・時刻・季節の生態学

**RSDSを生態学研究の一次保管先とする。** 旧ODSPから、公開原データに基づく太陽方位・日射・日周移動・季節的再配分・等日長比較・マクロ気候解析を引き継ぐ。一般的なODSP状態空間予測アルゴリズムは本研究の中心にしない。

## 現在の中心質問：時間利用と方向利用の組織化

> **太陽・季節・環境の変化に対して、生物の時間利用と移動方向はどのように組織化され、その違いは種や地域を越えて一般化するのか？**

日長・太陽方位・太陽高度・放射・気温・風を「環境状態」とし、観測された**昼夜の相対移動速度（時間ニッチのGPS代理量）**と、**太陽相対GPS移動方向（移動ニッチ）**を分ける。両者が同じ原調査で比較可能か、さらに別の原調査にも残るかを中心に検証する。撮影時刻は検出の観測過程を明示し、真の活動時間と区別する。

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
