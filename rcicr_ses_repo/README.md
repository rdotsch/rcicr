# 대상의 정보(SES)에 따른 정신적 표상 시각화: `rcicr` 기반 Reverse Correlation 파이프라인

이 저장소는 **noise-based reverse correlation**(2IFC) 방법을 사용해  
참가자들이 **낮은 SES(사회경제적 지위)**를 떠올릴 때 갖는 **얼굴의 정신적 표상(Classification Image; CI)**을 생성하고,
생성된 CI(및 anti-CI, z-map)를 재현 가능하게 정리한 프로젝트입니다.

- 핵심 방법: 무작위 노이즈가 덧씌워진(face **overlaid with random noise**) 얼굴 자극을 반복 제시하고(2IFC),
  참가자의 선택을 누적해 **정신적 표상의 시각적 프록시(CI)**를 산출합니다.
- 구현: Ron Dotsch의 **`rcicr` R 패키지**를 사용합니다.

> ⚠️ **중요(권한/윤리):**  
> 본 연구에서 사용한 **base face**(Moon et al., 2020 기반 한국인 합성 얼굴)는 연구자에게 직접 제공받은 자료입니다.  
> 따라서 이 저장소에는 **base face 원본 파일을 포함하지 않습니다.**  
> 재현을 위해서는 사용자가 적법한 경로로 base face를 확보해야 합니다.

---

## 1) 참고문헌 / 원 자료

- Brinkman, L., Todorov, A., & Dotsch, R. (2017). *Visualising mental representations: A primer on noise-based reverse correlation in social psychology.* **European Review of Social Psychology, 28**(1), 333–361. https://doi.org/10.1080/10463283.2017.1381469
- `rcicr` 패키지(원 저자 리포): https://github.com/rdotsch/rcicr
- `rcicr_examples`(예제 리포): https://github.com/rdotsch/rcicr_examples
- Moon et al. (2020) (base face 출처): 본 프로젝트에서는 **원본 base face를 배포하지 않음**(권한 제한)

---

## 2) 저장소 구조(권장)

```
rcicr_ses_repo/
├─ scripts/
│  ├─ 01_generate_stimuli.R
│  └─ 02_generate_classification_images.R
├─ data/
│  ├─ base_faces/          # (미포함) base face 이미지를 사용자가 넣는 폴더
│  └─ raw_responses/       # (미포함) Qualtrics/실험 로그 CSV를 넣는 폴더
├─ analysis/
│  ├─ input/               # (선택) 전처리된 CSV를 모으는 폴더
│  └─ output/              # (자동생성) CI 결과가 저장되는 폴더
├─ outputs/
│  ├─ ci/                  # 예시 CI/anti-CI 이미지
│  ├─ zmap/                # 예시 z-map 이미지
│  └─ poster/              # 포스터 PDF(선택)
└─ archive/                # 과거 버전 스크립트(정리용)
```

- `outputs/`에는 **최종 산출물 예시**(CI, z-map)를 포함할 수 있습니다.
- **원자료(참가자 응답 CSV)**는 개인정보/IRB 이슈가 있으므로, 공개 저장소에 올릴 때는 **익명화/요약본만** 권장합니다.

---

## 3) 설치(Installation)

### 3.1 R / RStudio
- R(권장: 최신 안정 버전)
- RStudio(선택이지만 권장)

### 3.2 패키지 설치
`rcicr`는 CRAN에서 제거된 이력이 있어(archived), **GitHub 설치를 권장**합니다.

```r
install.packages(c("remotes", "stringr", "dplyr"))
remotes::install_github("rdotsch/rcicr")
```

> Windows에서 설치 오류가 나면 Rtools가 필요할 수 있습니다.

---

## 4) 단계별 사용법(Quick start)

### Step 0. base face 준비
`data/base_faces/`에 base face 이미지를 넣습니다.

- 파일 조건(권장): 512×512, grayscale, 정면(face aligned)
- 예시 파일명(권장):  
  - `base_female.png`  
  - `base_male.png`

> ⚠️ 파일명은 스크립트에서 참조하므로, **스크립트 안의 base 파일명과 일치**시켜 주세요.

---

### Step 1. 자극(stimuli) 생성: `01_generate_stimuli.R`

이 스크립트는 2IFC 실험용 자극을 생성합니다. 핵심 함수는 `generateStimuli2IFC()`입니다.

```r
library(rcicr)

base <- list(
  "female" = "data/base_faces/base_female.png",
  "male"   = "data/base_faces/base_male.png"
)

generateStimuli2IFC(
  base_face_files = base,
  n_trials = 300,                 # 참가자 1인당 trial 수(예: 300)
  stimulus_path = "./stimuli",    # 생성 자극 저장 경로
  label = "ses",                  # 파일 prefix(실험에서 구분용)
  nscales = 5,
  noise_type = "gabor",
  sigma = 25
)
```

#### 주요 인자 설명(자주 바꾸는 것만)
- `n_trials`: trial 수(포스터/논문과 맞추기)
- `stimulus_path`: 자극 저장 폴더
- `label`: 파일명 prefix(실험 플랫폼에서 식별)
- `noise_type`: 일반적으로 `"gabor"`를 많이 씁니다
- `sigma`: 노이즈의 부드러움/스케일(문헌 및 파일럿 기반으로 조정)

실행 후 결과:
- `stimuli/` 폴더 아래에 이미지 세트(ori/inv)가 생성됨
- 동시에 `.Rdata` 파일이 생성됨(시드/노이즈 정보 포함)  
  → **CI 생성에 반드시 필요**하므로 보관해야 합니다.

---

### Step 2. 실험 실행 및 응답 수집(외부)
이 저장소는 “실험 플랫폼(예: Qualtrics/PsychoPy)” 자체를 포함하지 않습니다.  
대신, 실험을 통해 얻은 CSV 로그를 CI 생성 스크립트가 읽을 수 있도록 정리합니다.

#### CSV에 필요한 최소 컬럼(현재 스크립트 기준)
`02_generate_classification_images.R`는 아래 컬럼을 사용합니다.

- `id`: 참가자 ID
- `imgpath`: 자극 폴더 경로 문자열(예: `"stimuli_male/"`, `"stimuli_female/"`)
- `imgset`: 자극 세트 번호(= stim index)
- `selectedstim`: 참가자가 선택한 파일명(파일명에 `_ori` 또는 `_inv` 포함)
- (선택) `gender`, `age`, `nationality`, `ses.response`

> `selectedstim` 파일명 규칙: `..._ori.png` 또는 `..._inv.png`  
> 스크립트는 이를 이용해 선택을 +1/-1로 재코딩합니다.

---

### Step 3. CI 생성: `02_generate_classification_images.R`

이 스크립트는 여러 CSV를 모아서(참가자별/전체) CI를 생성합니다.
핵심 함수는 `generateCI2IFC()`입니다.

#### 실행 전 확인
- `analysis/`(또는 `data/raw_responses/`)에 CSV를 둡니다.
- Step 1에서 생성된 `.Rdata` 파일 경로를 스크립트에 지정합니다.
  - 예: `stimuli_female/rcic_seed_...Rdata`, `stimuli_male/rcic_seed_...Rdata`

#### 실행 후 산출물
- 전체 여성/남성 CI: `ci_all_F.png`, `ci_all_M.png` (예시 파일명)
- anti-CI(선택): `antici_all_F.png`, `antici_all_M.png`
- (선택) z-map 이미지: 픽셀 단위로 일관된 신호를 시각화

---

## 5) z-map(유의 픽셀 지도) 해석 가이드(중요)

z-map은 “어떤 픽셀이 선택과 일관되게 연결되는지”를 **시각화**해 줍니다.  
다만 픽셀 수가 매우 많기 때문에, z-map은 보통 다음처럼 취급합니다.

- ✅ **탐색적/시각적 보조 지표**로 사용(패턴 파악)
- ⚠️ 다중비교 문제 때문에, 이를 “확정적 특징”으로 단정하지 않기

현재 저장소의 예시 z-map:
- `outputs/zmap/ttest_base_image_ori1.png`
- `outputs/zmap/ttest_base_image_anti_ori1.png`

---

## 6) 최종 산출물(예시)

- `outputs/ci/`: 최종 CI / anti-CI 이미지 예시
- `outputs/zmap/`: z-map 예시
- `outputs/poster/`: 포스터 PDF(선택)

> 공개 GitHub에 올릴 때는, **참가자 원자료(raw data)**는 업로드하지 않는 것을 권장합니다.  
> 대신 CI 이미지(집계 결과)와 코드, 실행 방법을 공개하면 재현성/공유성이 크게 좋아집니다.

---

## 7) 자주 발생하는 오류 & 해결

### (1) `rcicr` 설치가 안 됨
- CRAN 설치가 아니라 GitHub 설치를 사용하세요:
  ```r
  remotes::install_github("rdotsch/rcicr")
  ```

### (2) 경로 문제(파일을 못 찾음)
- RStudio를 쓰면 스크립트 상단에서 `setwd()`를 명확히 지정하거나,
- 프로젝트(.Rproj)를 만들고 그 루트를 기준으로 경로를 통일하세요.

### (3) 한글 인코딩(cp949) 문제
- Windows/한글 환경에서 CSV 저장 시 `fileEncoding="cp949"`가 필요할 수 있습니다.
- Mac/Linux는 `UTF-8`로 저장 후 읽는 편이 안전합니다.

---

## 8) 인용(Citation)
본 저장소를 사용하거나 참고했다면 다음을 인용해 주세요.

- Brinkman, Todorov, & Dotsch (2017) primer (DOI 위 참고)
- `rcicr` GitHub repo: https://github.com/rdotsch/rcicr

---

## 9) 문의
- Maintainer: (여기에 본인 이름/메일)
