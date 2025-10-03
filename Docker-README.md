# Docker Setup for STP Challenge Test Script

This Docker setup is used to run the STP Challenge Baseline test script `test.py`.

## File Description

- **Dockerfile**: Container build configuration file
- **requirements.txt**: Python dependency packages list
- **docker-compose.yml**: Docker Compose configuration file (optional)

## Usage Instructions

### Method 1: Using Docker Command

#### 1. Build Image
```bash
docker build -t stp-challenge:latest .
```

#### 2. Run Container
```bash
# Ensure data directory exists
mkdir -p output

# Run test script
docker run --rm \
  -v $(pwd)/input:/app/input:ro \
  -v $(pwd)/output:/app/output \
  -v $(pwd)/checkpoints:/app/checkpoints:ro \
  stp-challenge:latest
```

### Method 2: Using Docker Compose (Recommended)

#### 1. Build and Run
```bash
# Automatically build and run
docker-compose up --build

# Run in background
docker-compose up -d --build
```

#### 2. View Logs
```bash
docker-compose logs -f
```

#### 3. Stop Container
```bash
docker-compose down
```

## Directory Structure

```
Project Root Directory/
├── input/           # Input data directory (read-only mount)
│   ├── valid_rna.h5ad
│   └── ...
├── output/          # Output results directory
│   └── submission.csv
├── checkpoints/     # Model checkpoint directory (read-only mount)
│   └── model_results.pkl
├── test.py         # Main test script
└── Dockerfile      # Docker configuration file
```

## Environment Requirements

- Docker 20.10+
- Docker Compose 1.29+ (optional)

## Notes

1. **Data Files**: Ensure `input/valid_rna.h5ad` and `checkpoints/model_results.pkl` files exist
2. **Output Directory**: Container will automatically create `output/` directory and generate `submission.csv`
3. **Permissions**: Ensure Docker has read permissions to the project directory
4. **Resources**: Recommend allocating at least 2GB memory to Docker

## Troubleshooting

### 1. Permission Issues
```bash
# Ensure directory permissions are correct
chmod -R 755 .
```

### 2. Insufficient Memory
```bash
# Check Docker resource usage
docker system df
```

### 3. Data Files Not Found
Ensure the following files exist:
- `input/valid_rna.h5ad`
- `checkpoints/model_results.pkl`

### 4. View Container Logs
```bash
# View logs from the last running container
docker logs $(docker ps -lq)

# View logs in real-time
docker logs -f <container_id>
```

## Development and Debugging

### Enter Container for Debugging
```bash
docker run -it --rm \
  -v $(pwd):/app \
  stp-challenge:latest \
  bash
```

### Rebuild (No Cache)
```bash
docker build --no-cache -t stp-challenge:latest .
```