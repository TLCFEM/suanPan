/*******************************************************************************
 * Copyright (C) 2017-2025 Theodore Chang
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

// ReSharper disable CppUseElementsView
#ifdef SUANPAN_VTK

// ReSharper disable StringLiteralTypo
#include "vtkParser.h"

#include <Domain/DomainBase.h>
#include <Element/Element.h>
#include <Toolbox/utility.h>
#include <vtkActor.h>
#include <vtkColorTransferFunction.h>
#include <vtkDataSetMapper.h>
#include <vtkDoubleArray.h>
#include <vtkInformation.h>
#include <vtkLookupTable.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkNamedColors.h>
#include <vtkPointData.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkScalarBarActor.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLMultiBlockDataWriter.h>

vtkInfo vtk_process(std::istringstream& command) {
    vtkInfo config;

    std::string keyword;

    while(!command.eof() && get_input(command, keyword))
        if(is_equal(keyword, "scale")) get_input(command, config.scale);
        else if(is_equal(keyword, "type") && get_input(command, keyword)) config.set(to_token(keyword));
        else if(is_equal(keyword, "fontsize")) get_input(command, config.font_size);
        else if(is_equal(keyword, "save")) get_input(command, config.file_name);
        else if(is_equal(keyword, "nobar")) config.color_bar = false;
        else if(is_equal(keyword, "element")) config.per_element = true;
        else if(is_equal(keyword, "material")) {
            config.per_material = true;
            config.per_section = false;
        }
        else if(is_equal(keyword, "section")) {
            config.per_section = true;
            config.per_material = false;
        }
        else if(is_equal(keyword, "size")) {
            get_input(command, config.canvas_size[0]);
            get_input(command, config.canvas_size[1]);
        }

    config.title_name = "Plotting Quantity " + std::string(to_name(config.display_type));

    return config;
}

void vtk_setup(vtkUnstructuredGrid* grid, const vtkInfo& config) {
    const vtkNew<vtkNamedColors> color;
    const vtkNew<vtkLookupTable> table;
    const vtkNew<vtkColorTransferFunction> func;
    const vtkNew<vtkDataSetMapper> mapper;
    const vtkNew<vtkScalarBarActor> bar;
    const vtkNew<vtkActor> actor;
    const vtkNew<vtkRenderWindowInteractor> interactor;
    const vtkNew<vtkRenderer> renderer;
    const vtkNew<vtkRenderWindow> window;

    func->SetColorSpaceToDiverging();
    func->AddRGBPoint(0., .230, .299, .754);
    func->AddRGBPoint(1., .706, .016, .150);
    table->SetNumberOfTableValues(256);
    table->Build();
    double rgb[4] = {0., 0., 0., 1.};
    for(auto I = 0; I < 256; ++I) {
        func->GetColor(static_cast<double>(I) / 256., rgb);
        table->SetTableValue(I, rgb);
    }

    const auto scalar = grid->GetPointData()->GetScalars();
    const auto index = std::max(0, std::min(to_index(config.display_type), scalar->GetNumberOfComponents() - 1));

    mapper->SetInputDataObject(grid);
    mapper->SelectColorArray(scalar->GetName());
    mapper->SetScalarModeToUsePointFieldData();
    mapper->SetArrayComponent(index);
    mapper->SetScalarRange(scalar->GetRange(index));
    mapper->SetLookupTable(table);

    bar->SetLookupTable(table);

    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(color->GetColor3d("DodgerBlue").GetData());
    actor->GetProperty()->EdgeVisibilityOn();
    actor->GetProperty()->SetOpacity(1.);
    actor->GetProperty()->SetLineWidth(4.);

    renderer->AddActor(actor);
    if(config.color_bar) renderer->AddViewProp(bar);
    renderer->SetBackground(color->GetColor3d("Grey").GetData());
    renderer->ResetCameraClippingRange();

    window->AddRenderer(renderer);

    interactor->SetRenderWindow(window);

    window->Render();
    window->SetSize(config.canvas_size[0], config.canvas_size[1]);
    window->SetWindowName(config.title_name.c_str());

    interactor->Start();
}

void vtk_save(vtkNew<vtkMultiBlockDataSet>&& grid, std::string file_name) {
    const vtkNew<vtkXMLMultiBlockDataWriter> writer;
    if(const auto ext = writer->GetDefaultFileExtension(); !file_name.ends_with(ext)) file_name += "." + std::string(ext);
    writer->SetInputData(grid);
    writer->SetFileName(file_name.c_str());
#ifdef SUANPAN_DEBUG
    writer->SetDataModeToAscii();
    writer->SetCompressorTypeToNone();
#else
    writer->SetDataModeToBinary();
#endif
    writer->Write();
    suanpan_debug("Plot is written to file \"{}\".\n", file_name);
}

int vtk_parser(const shared_ptr<DomainBase>& domain, std::istringstream& command) {
#ifdef SUANPAN_WIN
    domain->insert(std::async(std::launch::async, vtk_cell_plot, std::cref(domain), vtk_process(command)));
#else
    vtk_cell_plot(domain, vtk_process(command));
#endif

    return SUANPAN_SUCCESS;
}

class vtkBlock {
    vtkIdType num_node{0}, num_cell{0};

    const vtkNew<vtkDoubleArray> data;
    const vtkNew<vtkPoints> node;
    const vtkNew<vtkUnstructuredGrid> grid;

#ifdef SUANPAN_MT
    std::mutex mutex;
#endif

    mat tensor;
    u32_vec counter;

public:
    auto& init(const vtkIdType node_size, const vtkIdType cell_size, const char* name) {
        num_node = node_size;
        num_cell = cell_size;

        data->SetNumberOfComponents(6);
        data->SetNumberOfTuples(num_node);
        data->SetName(name);
        node->SetNumberOfPoints(num_node);

        for(auto I = 0; I < num_node; ++I) {
            node->SetPoint(I, 0., 0., 0.);
            data->SetTuple6(I, 0., 0., 0., 0., 0., 0.);
        }

        if(num_cell > 1) {
            tensor.zeros(6, num_node);
            counter.zeros(num_node);
        }

        grid->Allocate(num_cell);

        return *this;
    }

    auto add(const shared_ptr<Element>& element, const uvec& encoding, const vtkInfo& config) {
        const auto cell = element->Setup(encoding);
        if(!cell) return;

#ifdef SUANPAN_MT
        const std::scoped_lock lock(mutex);
#endif

        grid->InsertNextCell(cell->GetCellType(), cell->GetPointIds());

        const mat location = element->GetDeformation(config.scale).resize(3, encoding.n_elem);
        for(auto I = 0u; I < encoding.n_elem; ++I) node->SetPoint(static_cast<vtkIdType>(encoding(I)), location.colptr(I));

        auto result = element->GetData(config.record_type);
        if(result.empty()) return;

        result.resize(6, encoding.n_elem);
        if(num_cell > 1) {
            tensor.cols(encoding) += result;
            counter(encoding) += 1u;
        }
        else
            for(auto I = 0u; I < encoding.n_elem; ++I) data->SetTuple(static_cast<vtkIdType>(encoding(I)), result.colptr(I));
    }

    auto& attach() {
        if(num_cell > 1) {
            for(auto I = 0; I < num_node; ++I)
                if(counter(I) > 0u) {
                    tensor.col(I) /= counter(I);
                    data->SetTuple(I, tensor.colptr(I));
                }

            tensor.reset();
            counter.reset();
        }

        grid->SetPoints(node);
        grid->GetPointData()->SetScalars(data);

        return grid;
    }
};

void vtk_cell_single_block(const shared_ptr<DomainBase>& domain, vtkInfo&& config) {
    auto& node_map = domain->get_compact_node_map();
    auto& element_pool = domain->get_element_pool();

    vtkBlock block;
    block.init(static_cast<vtkIdType>(node_map.size()), static_cast<vtkIdType>(element_pool.size()), config.category.c_str());

    std::ranges::for_each(element_pool, [&](const shared_ptr<Element>& element) {
        auto encoding = element->get_node_encoding();
        for(auto& I : encoding) I = node_map.at(I);

        block.add(element, encoding, config);
    });

    if(config.file_name.empty()) vtk_setup(block.attach(), config);
    else {
        vtkNew<vtkMultiBlockDataSet> root;
        root->SetBlock(0u, block.attach());
        root->GetMetaData(0u)->Set(vtkCompositeDataSet::NAME(), "Default");
        domain->insert(std::async(std::launch::async, vtk_save, std::move(root), config.file_name));
    }
}

void vtk_cell_per_material(const shared_ptr<DomainBase>& domain, vtkInfo&& config) {
    if(config.file_name.empty()) return;

    std::unordered_map<uword, vtkIdType> element_count;

    auto& element_pool = domain->get_element_pool();
    for(auto& element : element_pool)
        for(const auto tag : element->get_material_tag()) element_count[tag] += 1;

    std::unordered_map<uword, vtkBlock> blocks;

    auto& compact_map = domain->get_compact_node_map_per_material();
    for(auto& [tag, node_map] : compact_map)
        if(element_count[tag] > 0) blocks[tag].init(static_cast<vtkIdType>(node_map.size()), element_count[tag], config.category.c_str());

    suanpan::for_all(element_pool, [&](const shared_ptr<Element>& element) {
        for(const auto tag : element->get_material_tag()) {
            auto encoding = element->get_node_encoding();
            auto& node_map = compact_map.at(tag);
            for(auto& I : encoding) I = node_map.at(I);

            blocks.at(tag).add(element, encoding, config);
        }
    });

    vtkNew<vtkMultiBlockDataSet> root;

    auto counter = 0u;
    for(auto& [tag, block] : blocks) {
        root->SetBlock(counter, block.attach());
        root->GetMetaData(counter++)->Set(vtkCompositeDataSet::NAME(), "Material " + std::to_string(tag));
    }

    domain->insert(std::async(std::launch::async, vtk_save, std::move(root), config.file_name));
}

void vtk_cell_per_section(const shared_ptr<DomainBase>& domain, vtkInfo&& config) {
    if(config.file_name.empty()) return;

    std::unordered_map<uword, vtkIdType> element_count;

    auto& element_pool = domain->get_element_pool();
    for(auto& element : element_pool)
        for(const auto tag : element->get_section_tag()) element_count[tag] += 1;

    std::unordered_map<uword, vtkBlock> blocks;

    auto& compact_map = domain->get_compact_node_map_per_section();
    for(auto& [tag, node_map] : compact_map)
        if(element_count[tag] > 0) blocks[tag].init(static_cast<vtkIdType>(node_map.size()), element_count[tag], config.category.c_str());

    suanpan::for_all(element_pool, [&](const shared_ptr<Element>& element) {
        for(const auto tag : element->get_section_tag()) {
            auto encoding = element->get_node_encoding();
            auto& node_map = compact_map.at(tag);
            for(auto& I : encoding) I = node_map.at(I);

            blocks.at(tag).add(element, encoding, config);
        }
    });

    vtkNew<vtkMultiBlockDataSet> root;

    auto counter = 0u;
    for(auto& [tag, block] : blocks) {
        root->SetBlock(counter, block.attach());
        root->GetMetaData(counter++)->Set(vtkCompositeDataSet::NAME(), "Section " + std::to_string(tag));
    }

    domain->insert(std::async(std::launch::async, vtk_save, std::move(root), config.file_name));
}

void vtl_cell_multiple_block(const shared_ptr<DomainBase>& domain, vtkInfo&& config) {
    if(config.file_name.empty()) return;

    suanpan::unordered_map<uword, vtkBlock> blocks;
    suanpan::for_all(domain->get_element_pool(), [&](const shared_ptr<Element>& element) {
        auto encoding = element->get_node_encoding();
        for(auto I = 0llu; I < encoding.n_elem; ++I) encoding(I) = I;

        blocks[element->get_tag()].init(static_cast<vtkIdType>(encoding.n_elem), 1, config.category.c_str()).add(element, encoding, config);
    });

    vtkNew<vtkMultiBlockDataSet> root;

    auto counter{0u};
    for(auto& [tag, block] : blocks) {
        root->SetBlock(counter, block.attach());
        root->GetMetaData(counter++)->Set(vtkCompositeDataSet::NAME(), "Element " + std::to_string(tag));
    }

    domain->insert(std::async(std::launch::async, vtk_save, std::move(root), config.file_name));
}

void vtk_cell_plot(const shared_ptr<DomainBase>& domain, vtkInfo config) {
    auto handler = vtk_cell_single_block;

    if(config.per_element) handler = vtl_cell_multiple_block;
    else if(config.per_material) handler = vtk_cell_per_material;
    else if(config.per_section) handler = vtk_cell_per_section;

    return handler(domain, std::move(config));
}

#endif
