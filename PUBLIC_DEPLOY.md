# Public Deployment

这个版本可以部署成公网 Web App，不依赖你的电脑一直开着。

## 推荐配置

RNA-seq / DESeq2 / clusterProfiler 比普通网页占内存多，建议云服务器至少：

- 2 CPU
- 4 GB RAM 起步，数据较大时建议 8 GB+
- Docker 支持

## 方案 A：Render

1. 把 `outputs/rna_seq_webapp` 目录上传到 GitHub 仓库。
2. 在 Render 创建 `New Web Service`。
3. 选择这个仓库。
4. Environment 选 `Docker`。
5. Dockerfile path 填：

```text
./Dockerfile
```

6. 等待 build 完成后，Render 会给一个 `https://...onrender.com` 地址。

如果使用 `render.yaml`，Render 也可以自动识别服务配置。

## 方案 B：Fly.io

在 `outputs/rna_seq_webapp` 目录中运行：

```bash
fly launch
fly deploy
```

如果 app 名字冲突，修改 `fly.toml` 里的 `app = "rnaseq-deseq-webapp"`。

## 方案 C：普通云服务器 / VPS

在服务器上安装 Docker，然后运行：

```bash
cd outputs/rna_seq_webapp
docker build -t rnaseq-deseq-webapp .
docker run -p 8788:8788 -e HOST=0.0.0.0 -e PORT=8788 rnaseq-deseq-webapp
```

然后访问：

```text
http://服务器IP:8788
```

如果要正式公开，建议再加 Nginx + HTTPS。

## 安全提醒

开放式应用意味着任何知道网址的人都可以上传文件并运行分析。正式公开前建议加：

- 登录密码
- HTTPS
- 上传文件大小限制
- 任务队列或并发限制
- 定期清理临时文件

