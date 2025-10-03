# STP Challenge Baseline - Docker Setup

Docker setup for running the STP Challenge Baseline test script `test.py`.

## Files

- **Dockerfile**: Container build configuration (Python 3.10)
- **requirements.txt**: Python dependencies with pinned versions
- **docker-compose.yml**: Docker Compose configuration

## Quick Start

### Method 1: Docker Compose (Recommended)

```bash
# Build and run
docker-compose up --build

# Run in background
docker-compose up -d --build

# View logs
docker-compose logs -f

# Stop
docker-compose down
```

### Method 2: Docker Command

#### 1. Build Image
```bash
docker build -t stp-baseline:latest .
```

#### 2. Run Container
```bash
# Ensure output directory exists
mkdir -p output

# Run test script
docker run --rm \
  -v $(pwd)/input:/app/input:ro \
  -v $(pwd)/output:/app/output \
  stp-baseline:latest
```

## Directory Structure

```
project_root/
├── input/                    # Required: Input data (mounted)
│   └── valid_rna.h5ad
├── output/                   # Generated: Results output
│   └── submission.csv
├── checkpoints/              # Built-in: Model files (included in image)
│   └── model_results.pkl
├── test.py                   # Main test script
├── train.py                  # Training script  
├── requirements.txt          # Python dependencies
├── Dockerfile               # Container config
└── docker-compose.yml       # Docker Compose config
```

## Prerequisites

- Docker 20.10+ 
- Docker Compose 1.29+ (optional)
- Input data: `input/valid_rna.h5ad`

## Configuration

The container includes:
- **Python 3.10** with scientific computing libraries
- **Pinned dependencies** matching your environment versions
- **Built-in models** - checkpoints embedded in image
- **Automatic output** generation to mounted directory

## Important Notes

1. **Input Required**: Only `input/valid_rna.h5ad` needed externally
2. **Models Included**: Checkpoints are built into the Docker image
3. **Output Generated**: Results saved to `./output/submission.csv`
4. **Memory**: Recommended 2GB+ for Docker

## Troubleshooting

### Permission Issues
```bash
# Fix directory permissions
chmod -R 755 .
```

### Missing Input Data
Ensure `input/valid_rna.h5ad` exists:
```bash
ls -la input/valid_rna.h5ad
```

### Memory Issues
```bash
# Check Docker resource usage
docker system df
docker stats
```

### View Logs
```bash
# Docker Compose
docker-compose logs -f

# Docker command
docker logs <container_id>
```

## Development

### Interactive Shell
```bash
# Enter container
docker run -it --rm \
  -v $(pwd)/input:/app/input:ro \
  -v $(pwd)/output:/app/output \
  stp-baseline:latest bash
```

### Clean Rebuild
```bash
# Remove old images and rebuild
docker-compose down
docker image rm stp-baseline:latest
docker-compose up --build --force-recreate
```