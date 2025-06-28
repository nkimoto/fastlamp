# Highly Optimized LCM makefile
# Enhanced version of: gcc -O3 -Os -s -o lcm lcm.c -lm

SAMPLE_SMALL = sample/sample_item.csv.4lcm53
SAMPLE_LARGE = sample/benchmark_item.csv.4lcm53

# コンパイラ設定
CC = gcc
MAKEFLAGS += -j$(shell nproc)

# ビルド最適化設定
TMPDIR ?= /tmp
export TMPDIR
ifneq ($(shell which ccache 2>/dev/null),)
    CC := ccache $(CC)
endif

# 基本最適化
CFLAGS_BASE = -O3 -march=native -mtune=native -pipe

# 数値計算最適化
CFLAGS_MATH = -ffast-math -fno-signed-zeros -fno-trapping-math -ffinite-math-only \
              -fno-rounding-math -fcx-limited-range

# ループ・ベクトル化最適化
CFLAGS_LOOP = -funroll-loops -fprefetch-loop-arrays -ftree-vectorize

# インライン・関数最適化
CFLAGS_INLINE = -finline-functions -fomit-frame-pointer \
                -ffunction-sections -fdata-sections

# アライメント最適化
CFLAGS_ALIGN = -falign-functions=32 -falign-loops=32 -falign-jumps=16

# 分岐・ブロック最適化
CFLAGS_BRANCH = -freorder-blocks

# メモリ最適化
CFLAGS_MEM = -fmerge-all-constants

# リンカ最適化
LDFLAGS = -s -Wl,-O2 -Wl,--gc-sections -Wl,--as-needed -Wl,--hash-style=gnu

# 統合フラグ
CFLAGS = $(CFLAGS_BASE) $(CFLAGS_MATH) $(CFLAGS_LOOP) $(CFLAGS_INLINE) \
         $(CFLAGS_ALIGN) $(CFLAGS_BRANCH) $(CFLAGS_MEM)

# 高速PGO最適化ビルド
lcm: lcm.c $(SAMPLE_SMALL) $(SAMPLE_LARGE)
	@echo "Building highly optimized LCM..."
	@# Phase 1: Profile generation
	@mkdir -p pgo
	$(CC) $(CFLAGS) -fprofile-generate=./pgo -o lcm_prof lcm.c -lm
	@# Phase 2: Intensive training for optimal profiling
	@echo "Running PGO training..."
	@./lcm_prof CIf $(SAMPLE_SMALL) 5 >/dev/null 2>&1 || true
	@./lcm_prof CIf $(SAMPLE_SMALL) 10 >/dev/null 2>&1 || true
	@./lcm_prof CIf $(SAMPLE_LARGE) 50 >/dev/null 2>&1 || true
	@./lcm_prof CIf $(SAMPLE_LARGE) 100 >/dev/null 2>&1 || true
	@./lcm_prof CIf $(SAMPLE_LARGE) 150 >/dev/null 2>&1 || true
	@# Phase 3: Final optimized build with collected profiles
	@echo "Building final optimized binary..."
	$(CC) $(CFLAGS) -fprofile-use=./pgo -fprofile-correction \
		$(LDFLAGS) -o $@ lcm.c -lm
	@rm -rf lcm_prof pgo
	@echo "✓ Highly optimized LCM built successfully!"

clean:
	@rm -rf lcm lcm_prof pgo *.gcda *.gcno
	@echo "✓ Cleaned up"

.PHONY: fast clean
