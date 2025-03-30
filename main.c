#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define CH 3
#define Ych 0
#define ROW 3
#define COL 3
#define IM_SIZE 1024

typedef struct {
    double x;
    double y;
    double r;
} Vec3d;

void get_data(void);
int cal_value(int lsb_index, int bytes);
void processing(void);
void put_data(void);
void rgb_to_ybr(void);
void ybr_to_rgb(void);

unsigned char header[54];
unsigned char*** imgin;
unsigned char*** imgin_rgb;
unsigned char*** imgout;
double GaussianFilter[IM_SIZE][IM_SIZE];

int width, height, alignment;

unsigned char*** allocateImgArray(void) {
    unsigned char*** array = (unsigned char***)malloc(CH * sizeof(unsigned char**));
    if (array == NULL) {
        perror("メモリ割り当てエラー (channels)");
        exit(1);
    }

    for (int i = 0; i < CH; i++) {
        array[i] = (unsigned char**)malloc(width * sizeof(unsigned char*));
        if (array[i] == NULL) {
            perror("メモリ割り当てエラー (width)");
            exit(1);
        }
        for (int j = 0; j < width; j++) {
            array[i][j] = (unsigned char*)malloc(height * sizeof(unsigned char));
            if (array[i][j] == NULL) {
                perror("メモリ割り当てエラー (height)");
                exit(1);
            }
        }
    }
    return array;
}

void freeImgArray(unsigned char*** array) {
    for (int i = 0; i < CH; i++) {
        for (int j = 0; j < width; j++) {
            free(array[i][j]);
        }
        free(array[i]);
    }
    free(array);
}

int main(void) {
    get_data();
    rgb_to_ybr();
    processing();
    ybr_to_rgb();
    put_data();
    return 0;
}

void get_data(void) {
    FILE* fp;
    char name[20];
    int filesize, offset, bits;
    int i, h, w;

    printf("入力ファイル名:");
    scanf("%s", name);
    fp = fopen(name, "rb");
    if (fp == NULL) {
        printf("\n%sをオープンできません.\n", name);
        exit(1);
    }
    else printf("\n%sをオープンしました.\n", name);

    fread(header, 1, 54, fp);

    filesize = cal_value(2, 4);
    offset = cal_value(10, 4);
    width = cal_value(18, 4);
    height = cal_value(22, 4);
    bits = cal_value(28, 2);
    alignment = filesize - offset - width * height * (bits / 8);

    imgin = allocateImgArray();
    imgin_rgb = allocateImgArray();
    imgout = allocateImgArray();

    for (h = height - 1; h >= 0; h--) {
        for (w = 0; w < width; w++) {
            for (i = 2; i >= 0; i--)
                imgin_rgb[i][w][h] = imgin[i][w][h] = fgetc(fp);
        }
    }

    fclose(fp);
    printf("%sをクローズしました.\n", name);
}

int cal_value(int lsb_index, int bytes) {
    int i;
    int value = header[lsb_index + bytes - 1];
    for (i = bytes - 2; i >= 0; i--) {
        value <<= 8;
        value += header[lsb_index + i];
    }
    return (value);
}

void createGaussianFilter(int blockSize, double sigma) {
    int n, m;
    int f_size = blockSize / 2;
    double sum = 0.0;

    for (n = -f_size; n <= f_size; n++) {
        for (m = -f_size; m <= f_size; m++) {
            GaussianFilter[m + f_size][n + f_size] = (1.0 / (2 * M_PI * sigma * sigma)) * exp(-((double)(m * m + n * n) / (2 * sigma * sigma)));
            sum += GaussianFilter[m + f_size][n + f_size];
        }
    }

    for (n = -f_size; n <= f_size; n++) {
        for (m = -f_size; m <= f_size; m++) {
            GaussianFilter[m + f_size][n + f_size] /= sum;
        }
    }
}

void applyGaussianFilter(unsigned char*** input, unsigned char*** output, int blockSize, double sigma) {
    int n, m, w, h, i, j;
    int f_size = blockSize / 2;
    double sum = 0.0;

    createGaussianFilter(blockSize, sigma);

    for (h = 0; h < height; h++) {
        for (w = 0; w < width; w++) {
            sum = 0.0;
            for (n = -f_size; n <= f_size; n++) {
                for (m = -f_size; m <= f_size; m++) {
                    i = w + m;
                    j = h + n;
                    if (i < 0) i = -i;
                    if (j < 0) j = -j;
                    if (i >= width) i = 2 * width - i - 1;
                    if (j >= height) j = 2 * height - j - 1;
                    sum += input[0][i][j] * GaussianFilter[m + f_size][n + f_size];
                }
            }
            output[0][w][h] = (unsigned char)(sum + 0.5);
            output[1][w][h] = input[1][w][h];
            output[2][w][h] = input[2][w][h];
        }
    }
}

void adaptiveThreshold(unsigned char*** input, unsigned char*** output, int blockSize, int C) {
    int h, w;
    unsigned char difference;

    applyGaussianFilter(input, output, blockSize, 3.0);

    for (h = 0; h < height; h++) {
        for (w = 0; w < width; w++) {
            difference = output[0][w][h] - C;
            if (imgin[0][w][h] > difference)
                output[0][w][h] = 255;
            else
                output[0][w][h] = 0;
            output[1][w][h] = input[1][w][h];
            output[2][w][h] = input[2][w][h];
        }
    }
}

// Helper function for comparing accumulator values
typedef struct { int* aux; } hough_cmp_gt_data;
hough_cmp_gt_data cmp_data;

int hough_cmp_gt(const void* a, const void* b) {
    int l1 = *(int*)a, l2 = *(int*)b;
    return cmp_data.aux[l1] > cmp_data.aux[l2] ? -1 : cmp_data.aux[l1] < cmp_data.aux[l2] ? 1 : l1 < l2 ? -1 : 1;
}

static void findLocalMaximums(int numrho, int numangle, int threshold, const int* accum, int* sort_buf, int* sort_buf_size) {
    *sort_buf_size = 0;
    for (int r = 0; r < numrho; r++) {
        for (int n = 0; n < numangle; n++) {
            int base = (n + 1) * (numrho + 2) + r + 1;
            if (accum[base] > threshold && accum[base] > accum[base - 1] && accum[base] >= accum[base + 1] &&
                accum[base] > accum[base - numrho - 2] && accum[base] >= accum[base + numrho + 2]) {
                sort_buf[(*sort_buf_size)++] = base;
            }
        }
    }
}

int Round(double value) { return (int)floor(value + 0.5); }

int HoughCircles(unsigned char*** input, Vec3d* detected_circles, int* numCircles, double dp, double minDist, int cannyThreshold, int accThreshold, int minRadius, int maxRadius) {
    *numCircles = 0;
    double idp = 1.0 / dp;

    short* dx = (short*)calloc(width * height, sizeof(short)), * dy = (short*)calloc(width * height, sizeof(short));
    unsigned char* edges = (unsigned char*)calloc(width * height, sizeof(unsigned char));
    if (!dx || !dy || !edges) {
        free(dx);
        free(dy);
        free(edges);
        fprintf(stderr, "Memory allocation failed.\n");
        return 0;
    }

    for (int y = 1; y < height - 1; y++) {
        for (int x = 1; x < width - 1; x++) {
            dx[y * width + x] = (short)((input[0][x + 1][y - 1] + 2 * input[0][x + 1][y] + input[0][x + 1][y + 1]) -
                (input[0][x - 1][y - 1] + 2 * input[0][x - 1][y] + input[0][x - 1][y + 1]));
            dy[y * width + x] = (short)((input[0][x - 1][y + 1] + 2 * input[0][x][y + 1] + input[0][x + 1][y + 1]) -
                (input[0][x - 1][y - 1] + 2 * input[0][x][y - 1] + input[0][x + 1][y - 1]));
        }
    }

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            unsigned char edge = (abs(dx[y * width + x]) + abs(dy[y * width + x]) > cannyThreshold) ? 255 : 0;
            edges[y * width + x] = edge;
            /*エッジ抽出画像出力用
            input[0][x][y] = edge;
            input[1][x][y] = 126;
            input[2][x][y] = 126;
            */
        }
    }

    int acols = Round(width * idp), arows = Round(height * idp), astep = acols + 2;
    int* accum = calloc((arows + 2) * (acols + 2), sizeof(int));
    if (!accum) {
        free(dx);
        free(dy);
        free(edges);
        fprintf(stderr, "Memory allocation failed.\n");
        return 0;
    }

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            if (edges[y * width + x]) {
                double vx = dx[y * width + x], vy = dy[y * width + x], mag = sqrt(vx * vx + vy * vy);
                if (mag < 1.0f) continue;
                int sx = Round((vx * idp) * 1024 / mag), sy = Round((vy * idp) * 1024 / mag);
                int x0 = Round(x * idp * 1024), y0 = Round(y * idp * 1024);
                for (int k1 = 0; k1 < 2; k1++, sx = -sx, sy = -sy) {
                    for (int x1 = x0 + minRadius * sx, y1 = y0 + minRadius * sy, r = minRadius;
                        r <= maxRadius; x1 += sx, y1 += sy, r++) {
                        int x2 = x1 >> 10, y2 = y1 >> 10;
                        if ((unsigned)x2 >= (unsigned)acols || (unsigned)y2 >= (unsigned)arows) break;
                        accum[y2 * astep + x2]++;
                    }
                }
            }
        }
    }

    int* centers_buf = malloc(sizeof(int) * (arows + 2) * (acols + 2));
    if (!centers_buf) {
        free(dx);
        free(dy);
        free(edges);
        free(accum);
        fprintf(stderr, "Memory allocation failed.\n");
        return 0;
    }

    int num_centers = 0;
    findLocalMaximums(acols, arows, accThreshold, accum, centers_buf, &num_centers);

    cmp_data.aux = accum;
    qsort(centers_buf, num_centers, sizeof(int), hough_cmp_gt);

    for (int i = 0; i < num_centers && *numCircles < 1024; i++) {
        int center_index = centers_buf[i], y = center_index / astep, x = center_index - y * astep;
        double curCenterX = (x + 0.5f) * dp, curCenterY = (y + 0.5f) * dp;
        double* dist_buf = malloc(sizeof(double) * width * height);
        if (!dist_buf) {
            free(dx);
            free(dy);
            free(edges);
            free(accum);
            free(centers_buf);
            fprintf(stderr, "Memory allocation failed.\n");
            return 0;
        }
        int dist_count = 0;
        for (int cy = 0; cy < height; cy++) {
            for (int cx = 0; cx < width; cx++) {
                if (edges[cy * width + cx]) {
                    double dx_c = curCenterX - cx, dy_c = curCenterY - cy, r2 = dx_c * dx_c + dy_c * dy_c;
                    if (r2 >= minRadius * minRadius && r2 <= maxRadius * maxRadius) {
                        dist_buf[dist_count++] = sqrt(r2);
                    }
                }
            }
        }

        if (dist_count > 0) {
            const int nBinsPerDr = 10, nBins = Round((maxRadius - minRadius) / dp * nBinsPerDr);
            int* bins = calloc(nBins, sizeof(int));
            if (!bins) {
                free(dist_buf);
                free(dx);
                free(dy);
                free(edges);
                free(accum);
                free(centers_buf);
                fprintf(stderr, "Memory allocation failed.\n");
                return 0;
            }
            for (int k = 0; k < dist_count; k++) {
                bins[max(0, min(nBins - 1, Round((dist_buf[k] - minRadius) / dp * nBinsPerDr)))]++;
            }

            int maxCount = 0;
            double best_r = 0;
            for (int j = nBins - 1; j > 0; j--) {
                if (bins[j]) {
                    int upbin = j, curCount = 0;
                    for (; j > upbin - nBinsPerDr && j >= 0; j--) {
                        curCount += bins[j];
                    }
                    double rCur = (upbin + j) / 2.0 / nBinsPerDr * dp + minRadius;
                    if ((curCount * best_r >= maxCount * rCur) || (curCount >= maxCount)) {
                        best_r = rCur;
                        maxCount = curCount;
                    }
                }
            }
            if (maxCount > accThreshold) {
                detected_circles[(*numCircles)++] = (Vec3d){ curCenterX, curCenterY, best_r };
            }
            free(bins);
        }
        free(dist_buf);
    }

    free(dx);
    free(dy);
    free(edges);
    free(accum);
    free(centers_buf);

    return 1;
}

// ユークリッド距離を計算する
double distance(Vec3d a, Vec3d b) {
    return sqrt(pow((a.x - b.x), 2) + pow((a.y - b.y), 2));
}

// DFSでグループ化する
void dfs(int index, int circle_count, Vec3d* circles, int* visited, int threshold, int* group_ids, int group_id) {
    visited[index] = 1;
    group_ids[index] = group_id;
    for (int i = 0; i < circle_count; i++) {
        if (!visited[i] && distance(circles[index], circles[i]) <= threshold) {
            dfs(i, circle_count, circles, visited, threshold, group_ids, group_id);
        }
    }
}

// グループごとの平均を計算して返す関数
Vec3d* compute_group_averages(Vec3d* circles, int circle_count, int threshold, int* group_count) {
    if (circle_count <= 0) {
        *group_count = 0;
        return NULL;
    }

    int* visited = (int*)calloc(circle_count, sizeof(int));
    int* group_ids = (int*)calloc(circle_count, sizeof(int));
    int current_group_id = 0;

    // グループを特定する
    for (int i = 0; i < circle_count; i++) {
        if (!visited[i]) {
            dfs(i, circle_count, circles, visited, threshold, group_ids, current_group_id);
            current_group_id++;
        }
    }

    *group_count = current_group_id;

    // 各グループの平均を計算
    Vec3d* group_averages = (Vec3d*)malloc(*group_count * sizeof(Vec3d));
    int* group_sizes = (int*)calloc(*group_count, sizeof(int));

    // 初期化
    for (int i = 0; i < *group_count; i++) {
        group_averages[i].x = 0.0f;
        group_averages[i].y = 0.0f;
        group_averages[i].r = 0.0f;
    }

    // グループごとの合計を計算
    for (int i = 0; i < circle_count; i++) {
        int group_id = group_ids[i];
        group_averages[group_id].x += circles[i].x;
        group_averages[group_id].y += circles[i].y;
        group_averages[group_id].r += circles[i].r;
        group_sizes[group_id]++;
    }

    // 合計から平均を計算
    for (int i = 0; i < *group_count; i++) {
        group_averages[i].x /= group_sizes[i];
        group_averages[i].y /= group_sizes[i];
        group_averages[i].r /= group_sizes[i];
    }

    // メモリ解放
    free(visited);
    free(group_ids);
    free(group_sizes);

    return group_averages;
}

void print_circle(unsigned char*** output, Vec3d circle, int thickness, int color_param) {
    double xy2, r2_inner, r2_outer;
    int w, h;

    r2_inner = pow(circle.r - thickness / 2, 2);
    r2_outer = pow(circle.r + thickness / 2, 2);

    for (h = (int)(circle.y - circle.r - thickness / 2 + 0.5); h < (int)(circle.y + circle.r + thickness / 2 + 0.5); h++) {
        for (w = (int)(circle.x - circle.r - thickness / 2 + 0.5); w < (int)(circle.x + circle.r + thickness / 2 + 0.5); w++) {
            xy2 = pow(w - circle.x, 2) + pow(h - circle.y, 2);
            if (xy2 >= r2_inner && xy2 <= r2_outer) {
                if (color_param == 1) {
                    output[0][w][h] = 145;
                    output[1][w][h] = 54;
                    output[2][w][h] = 34;
                }
                else if(color_param == 2){
                    output[0][w][h] = 170;
                    output[1][w][h] = 166;
                    output[2][w][h] = 16;
                }
                else {
                    output[0][w][h] = 231;
                    output[1][w][h] = 32;
                    output[2][w][h] = 143;
                }
            }
        }
    }
}

typedef struct {
    int x;
    int y;
} Point;

// キューの再利用用グローバル変数
Point* queue = NULL;
int queue_capacity = 0;

// 近傍の移動方向
const int dx[4] = { 1, -1, 0, 0 };
const int dy[4] = { 0, 0, 1, -1 };

void initialize_queue(int capacity) {
    if (queue != NULL) free(queue);
    queue = (Point*)malloc(capacity * sizeof(Point));
    queue_capacity = capacity;
}

int fill_th_area(unsigned char*** output, int start_x, int start_y, int t, int v) {
    int front = 0, rear = 0, cnt = 0;

    // 開始点をキューに追加
    queue[rear++] = (Point){ start_x, start_y };

    while (front < rear) {
        Point p = queue[front++];
        int w = p.x, h = p.y;

        // 範囲外チェックと値の確認
        if (w < 0 || w >= width || h < 0 || h >= height || output[0][w][h] != t) {
            continue;
        }

        // 塗りつぶし
        output[0][w][h] = v;
        cnt++;

        // 隣接ピクセルをキューに追加
        for (int i = 0; i < 4; i++) {
            int nw = w + dx[i];
            int nh = h + dy[i];
            queue[rear++] = (Point){ nw, nh };
        }
    }

    return cnt;
}

void molding(unsigned char*** output, int t) {
    initialize_queue(width * height);

    for (int h = 0; h < height; h++) {
        for (int w = 0; w < width; w++) {
            if (output[0][w][h] == 0) {
                if (fill_th_area(output, w, h, 0, 128) > t) {
                    fill_th_area(output, w, h, 128, 255);
                }
                else {
                    fill_th_area(output, w, h, 128, 0);
                }
            }
        }
    }

    free(queue);  // 最後に解放
}

void imgcpy(unsigned char*** input, unsigned char*** output) {
    int w, h;
    for (h = 0; h < height; h++) {
        for (w = 0; w < width; w++) {
            output[0][w][h] = input[0][w][h];
            output[1][w][h] = input[1][w][h];
            output[2][w][h] = input[2][w][h];
        }
    }
}

void get_ave_color(unsigned char*** input, Vec3d detected_circles, unsigned char* color_ave) {
    int w, h, count = 0;
    double r = 0, g = 0, b = 0;

    for (h = (int)(detected_circles.y - detected_circles.r + 0.5); h < (int)(detected_circles.y + detected_circles.r + 0.5); h++) {
        for (w = (int)(detected_circles.x - detected_circles.r + 0.5); w < (int)(detected_circles.x + detected_circles.r + 0.5); w++) {
            if (pow(w - detected_circles.x, 2) + pow(h - detected_circles.y, 2) < pow(detected_circles.r, 2)) {
                r += (double)input[0][w][h];
                g += (double)input[1][w][h];
                b += (double)input[2][w][h];
                count++;
            }
        }
    }

    color_ave[0] = (unsigned char)(r / count);
    color_ave[1] = (unsigned char)(g / count);
    color_ave[2] = (unsigned char)(b / count);
}

double get_pixel_sum(unsigned char*** input, Vec3d detected_circles) {
    int w, h;
    double count = 0;

    for (h = (int)(detected_circles.y - detected_circles.r + 0.5); h < (int)(detected_circles.y + detected_circles.r + 0.5); h++) {
        for (w = (int)(detected_circles.x - detected_circles.r + 0.5); w < (int)(detected_circles.x + detected_circles.r + 0.5); w++) {
            if (pow(w - detected_circles.x, 2) + pow(h - detected_circles.y, 2) <= pow(detected_circles.r, 2)) {
                if(input[0][w][h] == 0)
                    count++;
            }
        }
    }

    return count;
}

void RGBToHSV(double* RGB, double* HSV)
{
    double* H = &HSV[0], * S = &HSV[1], * V = &HSV[2];
    double R = RGB[0] / 255.0, G = RGB[1] / 255.0, B = RGB[2] / 255.0;

    // Calculate the max and min of red, green and blue.
    double mx = (R > G) ? ((R > B) ? R : B) : ((G > B) ? G : B);
    double mn = (R < G) ? ((R < B) ? R : B) : ((G < B) ? G : B);
    double delta = mx - mn;

    // Set the saturation and value
    *S = (mx != 0) ? (delta) / mx : 0;
    *V = mx;

    if (*S == 0.0f)
        *H = 0.0f;
    else
    {
        if (R == mx)        *H = (G - B) / delta;
        else if (G == mx)   *H = 2 + (B - R) / delta;
        else if (B == mx)   *H = 4 + (R - G) / delta;

        *H /= 6.0f;

        if (*H < 0.0f)      *H += 1.0f;
    }
}

void HSVToRGB(double* HSV, double* RGB)
{
    double H = HSV[0], S = HSV[1], V = HSV[2];
    double* R = &RGB[0], * G = &RGB[1], * B = &RGB[2];

    if (S == 0.0f)
    {
        // achromatic case
        *R = *G = *B = V * 255.0;
    }
    else
    {
        if (H == 1.0f)
            H = 0.0f;
        else
            H *= 6.0f;

        int i = (int)floor(H);
        double f = H - i;
        double p = V * (1 - S);
        double q = V * (1 - (S * f));
        double t = V * (1 - (S * (1 - f)));

        switch (i)
        {
        case 0: *R = V * 255.0; *G = t * 255.0; *B = p * 255.0; break;
        case 1: *R = q * 255.0; *G = V * 255.0; *B = p * 255.0; break;
        case 2: *R = p * 255.0; *G = V * 255.0; *B = t * 255.0; break;
        case 3: *R = p * 255.0; *G = q * 255.0; *B = V * 255.0; break;
        case 4: *R = t * 255.0; *G = p * 255.0; *B = V * 255.0; break;
        case 5: *R = V * 255.0; *G = p * 255.0; *B = q * 255.0; break;
        }
    }
}

int is_apple(unsigned char* input) {
    double hsv[3];
    double rgb[3] = {
        (double)input[0],
        (double)input[1],
        (double)input[2]
    };
    
    RGBToHSV(rgb, hsv);
    hsv[2] = 1.0;
    HSVToRGB(hsv, rgb);
    
    //オレンジ色の範囲内(下限、上限あり)なら柿orみかん。
    if (rgb[1] >= Round(rgb[2] * (50.0 / 71.0) + 90) && rgb[1] < Round(rgb[2] * (50.0 / 71.0) + 216))
        return 0;
    else
        return 1;
}

void processing(void) {
    unsigned char*** threshold = allocateImgArray();
    Vec3d detected_circles[512];
    Vec3d* group_averages;
    int numCircles = 0, group_count;
    double dp; int minDist, param1, param2, minRadius, maxRadius;
    unsigned char color_ave[CH];
    int apple = 0, orange = 0, persimmon = 0;

    imgcpy(imgin, imgout);
    adaptiveThreshold(imgin, threshold, 11, 3);
    molding(threshold, 11);

    
    if (HoughCircles(imgin, detected_circles, &numCircles,
        dp = 1.6, minDist = 100, param1 = 150, param2 = 100, minRadius = 50, maxRadius = 100)) {
        /*ハフ変換後の円描画用
        for (int i = 0; i < numCircles; i++)
            print_circle(imgout, detected_circles[i], 2, 1);
        */
        group_averages = compute_group_averages(detected_circles, numCircles, 100, &group_count);
        /*絞り込んだ円の描画用
        for (int i = 0; i < group_count; i++)
            print_circle(threshold, group_averages[i], 5, 1);
        imgcpy(threshold, imgout);
        */
        
        for (int i = 0; i < group_count; i++) {
            get_ave_color(imgin_rgb, group_averages[i], color_ave);
            double density = get_pixel_sum(threshold, group_averages[i]) / (M_PI * pow(group_averages[i].r, 2));
            if (!is_apple(color_ave)) {
                /*りんごの除外描画用
                print_circle(threshold, group_averages[i], 5, 1);
                */
                //double density = get_pixel_sum(threshold, group_averages[i]) / (M_PI * pow(group_averages[i].r, 2));
                if (density >= 0.1) {
                    orange++;
                    print_circle(imgout, group_averages[i], 7, 1);
                }
                else {
                    persimmon++;
                    print_circle(imgout, group_averages[i], 7, 2);
                }
            }
            else {
                apple++;
                print_circle(imgout, group_averages[i], 7, 0);
            }
        }
    }
    
    //デバッグ描画用
    //imgcpy(threshold, imgout);

    printf("\nりんご\t%d個\n", apple);
    printf("みかん\t%d個\n", orange);
    printf("柿\t%d個\n", persimmon);
    
    freeImgArray(threshold);
}

void rgb_to_ybr(void) {
    int i, h, w, ch, itemp;
    double dtemp[CH];
    double rgb_to_ybr[ROW][COL] = { { 0.2990, 0.5870, 0.1140},
                                   {-0.1687,-0.3313, 0.5000},
                                   { 0.5000,-0.4187,-0.0813} };

    for (h = 0; h < height; h++) {
        for (w = 0; w < width; w++) {
            for (ch = 0; ch < CH; ch++) {
                dtemp[ch] = 0.0;
                for (i = 0; i < COL; i++)
                    dtemp[ch] += rgb_to_ybr[ch][i] * (double)imgin[i][w][h];
            }

            for (ch = 0; ch < CH; ch++) {
                if (dtemp[ch] > 0.0)
                    itemp = (int)(dtemp[ch] + 0.5);
                else
                    itemp = (int)(dtemp[ch] - 0.5);
                if (ch != Ych)
                    itemp += 128;
                if (itemp > 255)
                    itemp = 255;
                else if (itemp < 0)
                    itemp = 0;
                imgin[ch][w][h] = itemp;
            }
        }
    }
}

void ybr_to_rgb(void) {
    int i, h, w, ch, itemp;
    double dtemp[CH];
    double ybr_to_rgb[ROW][COL] = { { 1.0000, 0.0000, 1.4020},
                                   { 1.0000,-0.3441,-0.7141},
                                   { 1.0000, 1.7720, 0.0000} };

    for (h = 0; h < height; h++) {
        for (w = 0; w < width; w++) {
            for (ch = 0; ch < CH; ch++) {
                dtemp[ch] = 0.0;
                for (i = 0; i < COL; i++) {
                    if (i == Ych)
                        dtemp[ch] += ybr_to_rgb[ch][i] * (double)imgout[i][w][h];
                    else
                        dtemp[ch] += ybr_to_rgb[ch][i] * (double)(imgout[i][w][h] - 128);
                }
            }

            for (ch = 0; ch < CH; ch++) {
                if (dtemp[ch] > 0.0)
                    itemp = (int)(dtemp[ch] + 0.5);
                else
                    itemp = (int)(dtemp[ch] - 0.5);
                if (itemp > 255)
                    itemp = 255;
                else if (itemp < 0)
                    itemp = 0;
                imgout[ch][w][h] = itemp;
            }
        }
    }
}

void put_data(void) {
    int i, w, h;
    FILE* fp;
    char name[20];

    printf("\n出力ファイル名:");
    scanf("%s", name);
    fp = fopen(name, "wb");
    if (fp == NULL) {
        printf("\n%sをオープンできません.\n", name);
        exit(1);
    }
    else printf("\n%sをオープンしました.\n", name);

    fwrite(header, 1, 54, fp);

    for (h = height - 1; h >= 0; h--) {
        for (w = 0; w < width; w++) {
            for (i = 2; i >= 0; i--)
                fputc(imgout[i][w][h], fp);
        }
    }

    for (i = 0; i < alignment; i++)
        fputc('\0', fp);

    freeImgArray(imgin);
    freeImgArray(imgin_rgb);
    freeImgArray(imgout);

    fclose(fp);
    printf("%sをクローズしました.\n", name);
}