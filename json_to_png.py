import json
import cv2
import numpy as np
import os

# 每张图像实例分割可视化
def save_instance_masks(img, anns, save_dir):
    mask = np.zeros((img['height'], img['width']), dtype=np.uint8)

    instance_id = 1
    for ann in anns:
        if ann['iscrowd'] == 1:
            continue

        pts = ann['segmentation'][0]
        pts = np.array(pts).reshape((-1, 2)).astype(np.int32)

        color = instance_id

        cv2.fillPoly(mask, [pts], color=color)
        instance_id += 1

    mask_file = os.path.join(save_dir, f"{img['file_name'][:-4]}_instance_masks.png")
    cv2.imwrite(mask_file, mask)


if __name__ == '__main__':
    json_file = 'LIVECell_dataset_2021/annotations/LIVECell/livecell_coco_test.json'
    save_dir = 'LIVECell_2021/test'

    if not os.path.exists(save_dir):
        os.mkdir(save_dir)

    with open(json_file, 'r') as f:
        data = json.load(f)

    for img in data['images']:
        anns = [ann for ann in data['annotations'] if ann['image_id'] == img['id']]
        save_instance_masks(img, anns, save_dir)
