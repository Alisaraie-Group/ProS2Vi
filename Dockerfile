FROM python:3.11-slim-bookworm

# Install system dependencies (DSSP, poppler, wkhtmltopdf, and rendering libraries)
RUN apt-get update && apt-get install -y --no-install-recommends \
    dssp \
    poppler-utils \
    wkhtmltopdf \
    # Font and rendering dependencies
    fontconfig \
    fonts-liberation \
    fonts-dejavu-core \
    libfontconfig1 \
    libfreetype6 \
    libx11-6 \
    libxext6 \
    libxrender1 \
    libxcb1 \
    libjpeg62-turbo \
    libpng16-16 \
    # SVG rendering support
    librsvg2-bin \
    libcairo2 \
    && rm -rf /var/lib/apt/lists/* \
    && fc-cache -f -v

# Set working directory
WORKDIR /app

# Install Python dependencies first (for caching)
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY src/ ./src/
COPY templates/ ./templates/
COPY pros2vi_gui.py .
COPY pros2vi_cli.py .

# Create output directory
RUN mkdir -p /app/output /app/uploads

# Set environment variable to skip browser auto-open
ENV DOCKER_ENV=1

# Expose Flask port
EXPOSE 3000

# Run the GUI application
CMD ["python", "pros2vi_gui.py"]
