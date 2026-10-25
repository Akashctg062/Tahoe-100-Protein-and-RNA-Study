# download_tier2.ps1 - Tier 2 open datasets for the extended (v2) study.  18 Oct 2026
# Run from PowerShell:
#   cd "E:\Protein and RNA\v2_extended\scripts"
#   powershell -ExecutionPolicy Bypass -File .\download_tier2.ps1            (steps 1-3, about 1.2 GB)
#   powershell -ExecutionPolicy Bypass -File .\download_tier2.ps1 -L1000     (adds step 4, LINCS L1000, several GB)
# Files go to  E:\Protein and RNA\data\raw\v2\  (new folder; existing data are not touched).
# curl.exe is built into Windows 10/11. "-C -" resumes an interrupted download: just run the script again.
param([switch]$L1000)
$ErrorActionPreference = "Stop"
$root = "E:\Protein and RNA\data\raw\v2"
New-Item -ItemType Directory -Force -Path $root, "$root\procan", "$root\ctrp", "$root\l1000", "$root\hansen2024" | Out-Null

function Get-File($url, $out) {
    Write-Host "`n>>> $out" -ForegroundColor Cyan
    curl.exe -L --fail --retry 5 --retry-delay 10 -C - -o $out $url
    if ($LASTEXITCODE -ne 0) { Write-Host "    FAILED (exit $LASTEXITCODE): $url" -ForegroundColor Red; return $false }
    $mb = [math]::Round((Get-Item $out).Length / 1MB, 1); Write-Host "    OK  $mb MB" -ForegroundColor Green; return $true
}

# 1) ProCan-DepMapSanger proteomics (Goncalves et al. 2022), Zenodo 6563157, ~592 MB
Get-File "https://zenodo.org/records/6563157/files/cancer_proteomics.zip?download=1" "$root\procan\cancer_proteomics.zip" | Out-Null

# 2) CTRP v2: try the original CTD2 release first, then the re-processed drevalpy copy (Zenodo 17600096, ~514 MB)
$ok = Get-File "https://ctd2-data.nci.nih.gov/Public/Broad/CTRPv2.0_2015_ctd2_ExpandedDataset/CTRPv2.0_2015_ctd2_ExpandedDataset.zip" "$root\ctrp\CTRPv2.0_2015_ctd2_ExpandedDataset.zip"
if (-not $ok) {
    Remove-Item "$root\ctrp\CTRPv2.0_2015_ctd2_ExpandedDataset.zip" -ErrorAction SilentlyContinue
    Get-File "https://zenodo.org/records/17600096/files/CTRPv2.zip?download=1" "$root\ctrp\drevalpy_CTRPv2.zip" | Out-Null
}

# 3) Small check list
Write-Host "`nSteps 1-2 finished. Now do the two manual/python steps described in the chat (ChEMBL hERG script; Hansen 2024 supplementary files)." -ForegroundColor Yellow

# 4) LINCS L1000 Level 5 (GEO GSE92742 + GSE70138) - only with -L1000; several GB
if ($L1000) {
    foreach ($gse in @(@("GSE92nnn","GSE92742"), @("GSE70nnn","GSE70138"))) {
        $base = "https://ftp.ncbi.nlm.nih.gov/geo/series/$($gse[0])/$($gse[1])/suppl/"
        Write-Host "`nListing $base" -ForegroundColor Cyan
        $html = (curl.exe -L -s $base) -join "`n"
        $names = [regex]::Matches($html, 'href="([^"]+)"') | ForEach-Object { $_.Groups[1].Value } |
                 Where-Object { $_ -match '(Level5_COMPZ.*\.gctx\.gz|sig_info.*\.txt\.gz|gene_info.*\.txt\.gz|cell_info.*\.txt\.gz|pert_info.*\.txt\.gz)$' } | Sort-Object -Unique
        if (-not $names) { Write-Host "    could not list files - please send me this message" -ForegroundColor Red }
        foreach ($n in $names) { Get-File ($base + $n) "$root\l1000\$n" | Out-Null }
    }
}
Write-Host "`nAll done. Folder: $root" -ForegroundColor Green
Get-ChildItem -Recurse $root | Where-Object { -not $_.PSIsContainer } | Select-Object FullName, @{n="MB";e={[math]::Round($_.Length/1MB,1)}} | Format-Table -AutoSize
