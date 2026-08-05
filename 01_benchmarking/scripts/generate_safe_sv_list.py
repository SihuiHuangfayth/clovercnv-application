import random
import sys

# --- 安全模式配置 (Safe Mode) ---
# 策略：降低密度，避开危险区，提高成功率

VARIANTS_PER_CHROM = 5        # 每条染色体只做 5 个
MIN_SIZE = 3000               
MAX_SIZE = 50000              # 最大 50kb
TYPES = ['DEL', 'DUP']        
BUFFER = 50000                # 两端留白 50kb
FAI_FILE = "/home/gaz6bs/data/row_Dd2/Dd2.fasta.fai"    

# ecDNA 配置
EC_DNA_RATIO = 0.3            
MIN_COPY = 5                  
MAX_COPY = 20                 

def parse_fai(fai_path):
    chroms = {}
    with open(fai_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            name = parts[0]
            length = int(parts[1])
            if "PfDd2" in name: 
                chroms[name] = length
    return chroms

def check_overlap(new_start, new_end, existing_regions):
    for start, end in existing_regions:
        if not (new_end < start or new_start > end):
            return True
    return False

def main():
    output_filename = sys.argv[1] if len(sys.argv) > 1 else "safe_sv_list.txt"

    try:
        chroms = parse_fai(FAI_FILE)
    except FileNotFoundError:
        print(f"Error: 找不到 {FAI_FILE}")
        sys.exit(1)

    print(f"正在生成安全版变异列表: {output_filename} ...")

    with open(output_filename, "w") as out:
        for chrom, length in chroms.items():
            if length < MAX_SIZE + 2 * BUFFER:
                continue

            existing_regions = []
            count = 0
            attempts = 0
            
            while count < VARIANTS_PER_CHROM and attempts < VARIANTS_PER_CHROM * 50:
                attempts += 1
                sv_type = random.choice(TYPES)
                size = random.randint(MIN_SIZE, MAX_SIZE)
                
                valid_max_start = length - BUFFER - size
                if valid_max_start <= BUFFER:
                    continue
                    
                start = random.randint(BUFFER, valid_max_start)
                end = start + size
                
                if check_overlap(start, end, existing_regions):
                    continue
                
                existing_regions.append((start, end))
                count += 1
                
                if sv_type == 'DUP':
                    if random.random() < EC_DNA_RATIO:
                        copy_num = random.randint(MIN_COPY, MAX_COPY)
                    else:
                        copy_num = 1
                    out.write(f"{chrom}\t{start}\t{end}\tDUP\t{copy_num}\n")
                else:
                    out.write(f"{chrom}\t{start}\t{end}\tDEL\n")

    print(f"成功! 已生成 {output_filename}")

if __name__ == "__main__":
    main()
