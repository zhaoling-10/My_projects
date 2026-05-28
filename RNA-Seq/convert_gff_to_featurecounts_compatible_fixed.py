"""
将 GFF3 格式转换为 featureCounts 兼容的 GTF 格式
保留原始的 gene_id 和 transcript_id
"""

import sys

def convert_gff_to_gtf(input_gff, output_gtf):
    # 第一遍扫描：建立 transcript -> gene 的映射关系
    transcript_to_gene = {}
    
    with open(input_gff, 'r') as infile:
        for line in infile:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            feature = fields[2]
            attributes = fields[8]
            
            # 解析属性
            attr_dict = {}
            for attr in attributes.split(';'):
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attr_dict[key] = value
            
            # 当遇到 mRNA/transcript 时，记录 transcript -> gene 关系
            if feature in ['mRNA', 'transcript']:
                transcript_id = attr_dict.get('ID', '')
                gene_id = attr_dict.get('Parent', '')
                if transcript_id and gene_id:
                    transcript_to_gene[transcript_id] = gene_id
    
    # 第二遍扫描：输出 GTF 格式
    with open(input_gff, 'r') as infile, open(output_gtf, 'w') as outfile:
        for line in infile:
            # 保留注释行
            if line.startswith('#'):
                outfile.write(line)
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            seqname, source, feature, start, end, score, strand, frame, attributes = fields
            
            # 解析 GFF 属性
            attr_dict = {}
            for attr in attributes.split(';'):
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attr_dict[key] = value
            
            # 跳过不需要的特征类型
            if feature not in ['gene', 'mRNA', 'transcript', 'exon', 'CDS', 'start_codon', 'stop_codon']:
                continue
            
            gene_id = None
            transcript_id = None
            
            # 处理 gene 特征
            if feature == 'gene':
                gene_id = attr_dict.get('ID', '')
                new_attributes = f'gene_id "{gene_id}"; gene_name "{gene_id}";'
            
            # 处理 mRNA/transcript 特征
            elif feature in ['mRNA', 'transcript']:
                transcript_id = attr_dict.get('ID', '')
                gene_id = attr_dict.get('Parent', '')
                new_attributes = f'gene_id "{gene_id}"; transcript_id "{transcript_id}";'
            
            # 处理 exon/CDS 等其他特征
            else:
                parent_id = attr_dict.get('Parent', '')
                
                # 从映射关系中查找 gene_id
                if parent_id in transcript_to_gene:
                    gene_id = transcript_to_gene[parent_id]
                    transcript_id = parent_id
                else:
                    # 如果找不到映射，尝试从 ID 推断
                    transcript_id = parent_id
                    # 例如：CM081674.1-g1.t1 -> CM081674.1-g1
                    if '.t' in parent_id:
                        gene_id = parent_id.rsplit('.t', 1)[0]
                    else:
                        gene_id = parent_id
                
                new_attributes = f'gene_id "{gene_id}"; transcript_id "{transcript_id}";'
            
            # 构建新的 GTF 行
            new_fields = [seqname, source, feature, start, end, score, strand, frame, new_attributes]
            outfile.write('\t'.join(new_fields) + '\n')

if __name__ == '__main__':
    if len(sys.argv) > 2:
        input_file = sys.argv[1]
        output_file = sys.argv[2]
    else:
        # 默认文件名
        input_file = 'Lepus_timidus_annotation.gff'
        output_file = 'genes_converted.gtf'
    
    print(f"开始转换 GFF 为 GTF...")
    print(f"输入文件: {input_file}")
    print(f"输出文件: {output_file}")
    convert_gff_to_gtf(input_file, output_file)
    print("✅ 转换完成！")
    print(f"验证转换结果：head -20 {output_file}")

