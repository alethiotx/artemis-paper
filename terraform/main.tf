terraform {
  # Backend configuration moved to backend.hcl
  # Initialize with: terraform init -backend-config=backend.hcl
  backend "s3" {}
}

# create ecr repo
resource "aws_ecr_repository" "this" {
  name = "kgpaper"
}

resource "aws_ecr_lifecycle_policy" "this" {
  repository = aws_ecr_repository.this.name

  policy = jsonencode(
    {
    "rules": [
      {
        "rulePriority": 1,
        "description": "Keep last 3 images",
        "selection": {
            "tagStatus": "any",
            "countType": "imageCountMoreThan",
            "countNumber": 3
        },
        "action": {
            "type": "expire"
        }
      }
    ]})
}

resource "aws_s3_bucket" "alethiotx_artemis" {
  bucket = "alethiotx-artemis"
}

# Account-level Public Access Block (must allow public bucket policies)
resource "aws_s3_account_public_access_block" "alethiotx_artemis" {
  block_public_acls       = true
  ignore_public_acls      = true
  block_public_policy     = false   # allow bucket policies that grant public access
  restrict_public_buckets = false
}

# Allow public reads of objects + optional list
resource "aws_s3_bucket_policy" "alethiotx_artemis" {
  bucket = aws_s3_bucket.alethiotx_artemis.id
  policy = jsonencode({
    Version = "2012-10-17",
    Statement = [
      {
        Sid       = "AllowPublicRead",
        Effect    = "Allow",
        Principal = "*",
        Action    = ["s3:GetObject"],
        Resource  = "${aws_s3_bucket.alethiotx_artemis.arn}/*"
      },
      {
        Sid       = "AllowPublicList",
        Effect    = "Allow",
        Principal = "*",
        Action    = ["s3:ListBucket"],
        Resource  = "${aws_s3_bucket.alethiotx_artemis.arn}"
      }
    ]
  })
  depends_on = [aws_s3_account_public_access_block.alethiotx_artemis]
}

# Bucket-level Public Access Block (must not block public policy)
resource "aws_s3_bucket_public_access_block" "alethiotx_artemis" {
  bucket                  = aws_s3_bucket.alethiotx_artemis.id
  block_public_acls       = true
  ignore_public_acls      = true
  block_public_policy     = false
  restrict_public_buckets = false
}

# CORS configuration for browser access
resource "aws_s3_bucket_cors_configuration" "alethiotx_artemis" {
  bucket = aws_s3_bucket.alethiotx_artemis.id

  cors_rule {
    allowed_headers = ["*"]
    allowed_methods = ["GET", "HEAD"]
    allowed_origins = ["*"]
    max_age_seconds = 3000
  }
}

# Website configuration
resource "aws_s3_bucket_website_configuration" "alethiotx_artemis" {
  bucket = aws_s3_bucket.alethiotx_artemis.id

  index_document {
    suffix = "index.html"
  }

  error_document {
    key = "error.html"
  }
}

resource "aws_s3_bucket" "alethiotx_artemis_internal" {
  bucket = "alethiotx-artemis-internal"
}

# Public ECR repository (ECR Public)
resource "aws_ecrpublic_repository" "this" {
  provider = aws.us_east_1
  repository_name = "artemis-paper"

  catalog_data {
    about_text = "Public image for Artermis paper repository."
  }
}

# Allow anonymous pull
resource "aws_ecrpublic_repository_policy" "this" {
  provider        = aws.us_east_1
  repository_name = aws_ecrpublic_repository.this.repository_name

  policy = jsonencode({
    Version = "2012-10-17",
    Statement = [{
      Sid       = "AllowPublicPull"
      Effect    = "Allow"
      Principal = "*"
      Action = [
        "ecr-public:BatchCheckLayerAvailability",
        "ecr-public:BatchGetImage",
        "ecr-public:GetDownloadUrlForLayer"
      ]
    }]
  })
}

output "public_image_uri_latest" {
  value = "${aws_ecrpublic_repository.this.repository_uri}:latest"
}