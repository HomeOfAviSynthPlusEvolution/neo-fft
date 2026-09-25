# neo-fft

[English](README.md) | [简体中文](README.zh-CN.md) | **日本語**

neo-fft は VapourSynth と AviSynth 用の周波数領域フィルタープラグインです。Neo FFT3D と Neo DFTTest を独立に再実装して一つのプラグインにまとめ、空間・時間方向のノイズ除去、Kalman フィルタリング、周波数領域でのシャープ化、ハロー除去、調整可能なスペクトルフィルタリングを提供します。

実装は C++17 を使用し、スカラーカーネル、Google Highway によるクロスプラットフォーム SIMD、PocketFFT に基づく変換を備えています。DualSynth2 が計算コアを両ホストに接続します。VapourSynth では `core.neo_fft`、AviSynth では `neo_fft_` 接頭辞の関数を使用します。

## 設計

neo-fft は周波数領域の計算とホストのフレーム管理を分離しています。コアはブロックの取り出し、窓掛け、FFT、スペクトルフィルタリング、再構成を担当し、ホスト層はパラメーター、フレーム要求、プロパティ、出力の割り当てを管理します。FFT3D と DFTTest は基盤を共有しつつ、それぞれの窓、ノイズモデル、オーバーラップ加算の規則を維持します。コアは単独でビルドおよびテストできます。

実装は動作仕様に基づいて開発されています。スカラー実装と独立した数学的定義で計算を確認し、SIMD 経路および固定した参照プラグインとの公開動作の比較を検証します。同じパラメーターでも、過去のすべての FFT3DFilter、DFTTest、Neo ビルドと画素単位で同一の出力を保証するものではありません。浮動小数点の丸め、しきい値による分岐、Kalman の状態管理、ディザリングは結果に影響する場合があります。具体的な差異は[移行ガイド（英語）](docs/api/en/migration.md)を参照してください。

両フィルターと FFT はホストの呼び出しスレッド上で処理を行い、ワーカースレッドやスレッドプールを作成しません。ホストは複数フレームを並行して要求できます。SIMD やバッチ FFT は内部のマルチスレッド処理を意味しません。

## 対応する処理

| 関数 | 用途 |
|---|---|
| `FFT3D` | 単一フレームまたは 2–5 フレームの Wiener ノイズ除去、Kalman フィルタリング、シャープ化、ハロー除去。周波数依存およびサンプリングによるノイズモデルに対応。 |
| `DFTTest` | 5 種類のスペクトルフィルター、空間・時間方向のオーバーラップ加算、12 種類の窓、周波数曲線、ノイズサンプリング、8 ビット出力のディザリング。 |
| `KernelInfo` | 自動選択された SIMD ターゲット、FFT 実装、ベクトル幅の取得。 |

両映像フィルターは、形式とサイズが固定されたプレーナー GRAY/YUV/RGB の 8/10/12/14/16 ビット整数および 32 ビット浮動小数点サンプルに対応します。AviSynth ではプレーナー YUVA/RGBA にも対応します。出力のサイズ、形式、フレーム数、フレームレートは入力と同じです。ブロック形状と境界の条件は、処理する各プレーンで確認されます。

既定ではアルファ以外の全プレーンを処理し、アルファはコピーします。`planes=[0]` は最初のプレーンだけを処理し、`planes=[]` はどのプレーンも処理しません。AviSynth は従来の `y/u/v/a` プレーンモードも受け付けます。モード 1 は出力プレーンに書き込まないことを明示し、後でそのプレーンを破棄するスクリプト向けです。優先順位は API リファレンスを参照してください。

FFT3D の既定値は `bt=3` で、隣接する 3 フレームを使用します。`bt=1` は空間方向のみのノイズ除去、`bt=0` は Kalman 処理です。Kalman は既定で最大 8 枚の過去フレームをウォームアップに使用し、近くのチェックポイントも再利用するため、シークのたびに映像の先頭から再計算しません。キャッシュの内容と要求履歴は再帰計算の結果に影響する場合があります。通常の時間フィルタリング用の未処理スペクトルキャッシュは、既定の予算が 128 MiB で、`cache_mb` と `cache_frames` で調整できます。これはプロセス全体のメモリ上限ではありません。

DFTTest の既定値は `tbsize=1` で、現在のフレームだけを使用します。`tbsize` を増やすと時間方向のフィルタリングが有効になり、`tmode` で中心出力または時間方向のオーバーラップ加算を選択します。どちらのフィルターも動きベクトルの推定や動き補償は行いません。

## ドキュメントと使用方法

API リファレンスでは関数の呼び出し方を、ナレッジベースでは窓、変換、スペクトルモデル、再構成によって入力から出力を求める過程を解説します。英語版を参照してください。

- [API reference (English)](docs/api/en/README.md)：関数シグネチャ、パラメーター、既定値、使用例。
- [Knowledge base (English)](docs/knowledge/en/README.md)：データ表現、数式、演算順序、境界、精度。
- [移行ガイド（英語）](docs/api/en/migration.md)：旧 Neo 系フィルターから移行する際に必要な変更と出力の差異。

ビルドしたプラグインを明示的に読み込むか、VapourSynth のプラグイン自動読み込みディレクトリに配置してください。以下は Windows のファイル名を使用しています。Linux では `neo-fft.so` を使用し、他のプラットフォームでも実際のプラグインパスに置き換えてください。VapourSynth のプラグイン識別子は `org.neofilters.neo_fft` です。

```python
import vapoursynth as vs

core = vs.core
core.std.LoadPlugin(path="/path/to/neo-fft.dll")

print(core.neo_fft.KernelInfo())
clip = core.std.BlankClip(width=640, height=360, format=vs.YUV420P8, length=24)
output = core.neo_fft.FFT3D(clip, sigma=2.0, bt=3, planes=[0])
# 別の選択肢：元の入力に DFTTest の空間ノイズ除去を適用。
# output = core.neo_fft.DFTTest(clip, sigma=8.0, tbsize=1, planes=[0])
output.set_output()
```

この最小例は合成クリップで呼び出し方を示します。両フィルターでは `sigma` の数学的な意味が異なるため、同じ値を指定してもノイズ除去の強さは一致しません。それぞれ個別に調整してください。

同じプラグインファイルに AviSynth C++ インターフェースも含まれます。インターフェースのバージョン 11 に対応したホストで `LoadPlugin` を使用してください。

```avs
LoadPlugin("/path/to/neo-fft.dll")
clip = BlankClip(width=640, height=360, length=24, pixel_type="YV12")
return neo_fft_FFT3D(clip, sigma=2.0, bt=3, planes=[0])
# DFTTest を使う場合：
# return neo_fft_DFTTest(clip, sigma=8.0, tbsize=1, planes=[0])
```

AviSynth の配列パラメーターは `[0, 1]` などのネイティブ配列を受け取り、単一の値も 1 要素の配列として扱います。DFTTest の曲線とノイズサンプリング位置は、空白、カンマ、コロンで区切った数値文字列も受け付けます。両フィルターは入力の音声とフィールドパリティを引き継ぎ、出力フレームのプロパティは対応する入力フレームから取得します。

パラメーターの省略を明示するには、VapourSynth では `None`、AviSynth では `Undefined()` を使用します。たとえば `planes=None` または `planes=Undefined()` は既定のプレーン選択を使用し、後者では AVS の `y/u/v/a` も有効になります。空配列 `[]` は処理対象のプレーンがないことを明示します。

FFT3D の `mt/ncpu/measure`、DFTTest の `threads/fft_threads`、両フィルターの `fft_backend` は互換性のために受け付けますが、内容は完全に無視します。スレッドを作成したり FFT バックエンドを選択したりすることはありません。旧スクリプトの移行では名前付き引数を推奨します。パラメーター順と移行時の制約は各 API を参照してください。

## SIMD と CPU 選択

SIMD を有効にしたビルドは、実行中の CPU が対応し、かつビルドに含まれる Highway ターゲットを選択します。スカラーへのフォールバックも利用できます。FFT は独立に実装を選択し、一部の固定サイズには専用の変換経路もあります。

両フィルターの `opt=1` は本プロジェクトのスカラーカーネルを選択しますが、FFT をスカラーに固定するものではありません。他の対応値では自動 SIMD を使用し、従来の `opt` の値で ISA を指定・制限することはできません。ビルド全体で SIMD を無効にするには `NEO_FFT_ENABLE_SIMD=OFF` を指定してください。

`core.neo_fft.KernelInfo()` は `fft_backend`、`target`、`fft`、`fft_lanes` を含む辞書を返し、フィールド名で参照します。AviSynth の `neo_fft_KernelInfo()` は `[fft_backend, target, fft, fft_lanes]` の固定順序で配列を返します。この問い合わせは自動選択の情報であり、特定のフィルターインスタンスの実行履歴ではありません。`fft_lanes` はベクトルのレーン数であり、スレッド数や高速化率ではありません。

FFT ターゲットは一般カーネルのターゲットと異なる場合があります。広い SIMD が常に高いスループットを保証するわけではありません。[KernelInfo（英語）](docs/knowledge/en/kernel-info.md)と[実行と精度（英語）](docs/knowledge/en/shared/execution-precision.md)を参照してください。

## ビルドとテスト

CMake 3.24 以降、Git、C++17 対応コンパイラーが必要です。CMake は固定バージョンの DualSynth2 と PocketFFT を取得し、SIMD 有効時には Highway 1.4.0 も取得します。両ホストの SDK はローカルで検出するか、自動取得します。

以下の手順は、映像ホストをインストールせずに両ホスト用プラグインとコアテストをビルドします。Windows では既定の AVS C++ インターフェースに MSVC または clang-cl が必要です。MinGW ビルドでは AVS インターフェースを無効にしてください。

```sh
cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DNEO_FFT_TEST_VAPOURSYNTH=OFF
cmake --build build/release --config Release --parallel 4
ctest --test-dir build/release -C Release --output-on-failure
```

| オプション | 用途 |
|---|---|
| `NEO_FFT_BUILD_VAPOURSYNTH=OFF` | VapourSynth のエントリーポイントを無効化。 |
| `NEO_FFT_BUILD_AVISYNTH=OFF` | AviSynth のエントリーポイントを無効化。両ホストを OFF にするとコアのみをビルド。 |
| `NEO_FFT_ENABLE_SIMD=OFF` | SIMD カーネルとベクトル化 FFT を無効化。 |
| `BUILD_TESTING=OFF` | テストをビルドしない。 |
| `NEO_FFT_TEST_VAPOURSYNTH=OFF` | コアテストとプラグインを残し、VapourSynth ホストテストを省略。VS インターフェースとテストが有効な場合、既定は ON。 |
| `Python3_EXECUTABLE=/path/to/python` | ホストテスト用 Python を指定。VS テストには VapourSynth をインポートでき、同じアーキテクチャのランタイムを読み込める環境が必要。 |
| `NEO_FFT_VS_SDK=/path/to/sdk` | ローカルの VapourSynth SDK を指定。 |
| `NEO_FFT_AVS_SDK=/path/to/sdk` | ローカルの AviSynth SDK を指定。 |
| `NEO_FFT_TEST_AVISYNTH=ON` | AviSynth ホストテストを有効化。既定では無効。Python と同じアーキテクチャのランタイムが必要。 |
| `NEO_FFT_AVISYNTH_RUNTIME=/path/to/avisynth.dll` | AviSynth ホストテスト用のランタイムライブラリを指定。 |
| `FETCHCONTENT_SOURCE_DIR_DUALSYNTH2=/path/to/dualsynth2` | 固定バージョンの取得に代えてローカルの DualSynth2 ソースを使用。 |

プラグインターゲット名は `neo_fft`、出力ファイルの基本名は `neo-fft` で、既定では両ホストのエントリーポイントを含みます。テストは FFT、窓、フィルター演算、境界、スカラー/SIMD 比較、ホスト動作を扱います。旧プラグインとのブラックボックス比較には、固定した参照バイナリーが別途必要です。

CI は Windows x64、Linux x64、macOS ARM64、Linux ASan/UBSan の検査を設定しています。Linux runner は Ubuntu 26.04 とシステム既定の GCC を使用し、Sanitizer 検査には Clang 22 を使用します。リリースワークフローは Windows、Linux、macOS の x64/ARM64 産物をビルドします。VapourSynth ホストテストは現在 Windows x64 で実行し、AviSynth ホストテストは上記のオプションで別途有効にします。パッケージには実際のビルド環境とテスト範囲を記録しており、すべての Linux ディストリビューションでの互換性を示すものではありません。

## 性能

既存の計測では、FFT3D の一般的な処理経路で旧 Neo FFT3D の約 **1.72–3.05 倍**、DFTTest の一般的な Wiener 経路で旧 Neo DFTTest の約 **4.58–7.04 倍**のスループットを得ています。比率は **neo-fft のスループット / 参照フィルターのスループット**、すなわち参照フィルターの時間 / neo-fft の時間です。**1 より大きい場合は neo-fft が高速**です。以下の範囲は 8/16 ビット整数、32 ビット浮動小数点、AVX2/AVX-512 の各設定をまとめたもので、信頼区間ではありません。

| 主な処理経路 | 相対スループット |
|---|---:|
| FFT3D 空間ノイズ除去（`bt=1`） | 2.39–2.85× |
| FFT3D 2 フレームノイズ除去（`bt=2`） | 2.32–3.05× |
| FFT3D 3 フレームノイズ除去（`bt=3`、既定の時間モード） | 1.97–2.53× |
| FFT3D 4 フレームノイズ除去（`bt=4`） | 1.86–2.37× |
| FFT3D 5 フレームノイズ除去（`bt=5`） | 1.72–2.15× |
| FFT3D Kalman（`bt=0`、順次要求） | 1.90–2.43× |
| DFTTest 空間 Wiener（`tbsize=1`） | 4.60–6.65× |
| DFTTest 3 フレーム Wiener（`tbsize=3, tmode=0`） | 5.25–7.04× |
| DFTTest 5 フレーム Wiener（`tbsize=5, tmode=0`） | 4.58–5.77× |

この範囲は、Highway で近代化された旧 Neo FFT3D / Neo DFTTest を参照とする、AviSynth 上の単一スレッドでの過去の比較結果です。以降の最適化は計測に含まれず、処理チェーン全体のスループットを示すものでもありません。実際の結果は入力、パラメーター、ハードウェア、ホストの並行処理数によって変わります。

## 開発と貢献

メンテナーが技術方針、変更のレビュー、リリースを担当します。不具合報告、提案、貢献を歓迎します。数値動作、公開インターフェース、重要な設計変更については、実装前に目的と方針を相談してください。

本プロジェクトでは実装、テスト、レビューに AI を活用します。貢献には問題、方法、検証内容、AI の関与を記載してください。報告にはバージョン、OS、CPU、コンパイラー、ビルド設定、入出力形式、最小再現例を含めてください。数値差分の報告には参照バージョン、パラメーター、要求順序を、性能報告には計測範囲とスレッド設定も記載してください。

## 謝辞とライセンス

neo-fft のインターフェースと周波数領域のフィルタリング機能の基礎を築いた、以下の上流プロジェクトの作者と貢献者の皆様に感謝します。

- [FFT3DFilter](https://github.com/pinterf/fft3dfilter)：Alexander G. Balakhnin（Fizick）が開発し、martin53 が AviSynth 2.6 に対応させ、Ferenc Pintér（pinterf）が高ビット深度対応などの改良を加えました。
- [DFTTest](https://github.com/pinterf/dfttest)：tritical が開発し、Firesledge が 16 ビット処理を追加、DJATOM が AviSynth+ に移植し、pinterf が高ビット深度対応とクロスプラットフォーム対応を進めました。

neo-fft は以下のライブラリも使用しています。

- [Google Highway](https://github.com/google/highway)：クロスプラットフォームの SIMD を提供します。
- [PocketFFT](https://github.com/mreineck/pocketfft)：FFT3D と DFTTest の周波数領域変換に使用します。
- [DualSynth2](https://github.com/HomeOfAviSynthPlusEvolution/dualsynth2)：VapourSynth、AviSynth と共有計算コアを接続します。

テスト、問題報告、改善に協力する開発者とユーザーの皆様に感謝します。

開発に使用する LLM サブスクリプションをご支援いただいた [SB.SB](https://sb.sb) に感謝します。

neo-fft は GNU General Public License バージョン 2 以降（`GPL-2.0-or-later`）で提供します。全文は [LICENSE](LICENSE) を参照してください。第三者コンポーネントはそれぞれの著作権表示とライセンス条項を維持します。
