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

#ifdef SUANPAN_VTK

// ReSharper disable StringLiteralTypo
#include "vtkParser.h"

#include <Domain/DomainBase.h>
#include <Domain/Node.h>
#include <Element/Element.h>
#include <Toolbox/utility.h>
#include <vtkAutoInit.h>
VTK_MODULE_INIT(vtkRenderingOpenGL2)  // NOLINT(cppcoreguidelines-special-member-functions, hicpp-special-member-functions)
VTK_MODULE_INIT(vtkInteractionStyle)  // NOLINT(cppcoreguidelines-special-member-functions, hicpp-special-member-functions)
VTK_MODULE_INIT(vtkRenderingFreeType) // NOLINT(cppcoreguidelines-special-member-functions, hicpp-special-member-functions)

#include <vtkActor.h>
#include <vtkColorTransferFunction.h>
#include <vtkDataSetMapper.h>
#include <vtkDoubleArray.h>
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
#include <vtkXMLUnstructuredGridWriter.h>

vtkInfo vtk_process(std::istringstream& command) {
    vtkInfo config;

    std::string keyword;

    while(!command.eof() && get_input(command, keyword))
        if(is_equal(keyword, "scale") && !get_input(command, config.scale)) config.scale = 1.;
        else if(is_equal(keyword, "deformed")) config.on_deformed = true;
        else if(is_equal(keyword, "undeformed")) config.on_deformed = false;
        else if(is_equal(keyword, "type") && get_input(command, keyword)) config.type = to_token(keyword);
        else if(is_equal(keyword, "fontsize") && !get_input(command, config.font_size)) config.font_size = 8;
        else if(is_equal(keyword, "save") && get_input(command, config.file_name)) config.save_file = true;
        else if(is_equal(keyword, "nobar")) config.colorbar = false;
        else if(is_equal(keyword, "noaverage")) config.average = false;
        else if(is_equal(keyword, "material") && !get_input(command, config.material_type)) config.material_type = -1;
        else if(is_equal(keyword, "size")) {
            if(!get_input(command, config.canvas_size[0])) config.canvas_size[0] = 500;
            if(!get_input(command, config.canvas_size[1])) config.canvas_size[1] = 500;
        }

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

    mapper->SetInputDataObject(grid);
    mapper->SetLookupTable(table);
    mapper->SetScalarRange(grid->GetPointData()->GetScalars()->GetRange());

    bar->SetLookupTable(mapper->GetLookupTable());

    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(color->GetColor3d("DodgerBlue").GetData());
    actor->GetProperty()->EdgeVisibilityOn();
    actor->GetProperty()->SetOpacity(1.);

    renderer->AddActor(actor);
    if(config.colorbar) renderer->AddViewProp(bar);
    renderer->SetBackground(color->GetColor3d("Grey").GetData());
    renderer->ResetCameraClippingRange();

    window->AddRenderer(renderer);

    interactor->SetRenderWindow(window);

    window->Render();
    window->SetSize(config.canvas_size[0], config.canvas_size[1]);
    window->SetWindowName(config.title_name.c_str());

    interactor->Start();
}

void vtk_save_single(vtkNew<vtkUnstructuredGrid>&& grid, std::string file_name) {
    // NOLINT(performance-unnecessary-value-param)
    const vtkNew<vtkXMLUnstructuredGridWriter> writer;
    if(const auto ext = writer->GetDefaultFileExtension(); !file_name.ends_with(ext)) file_name += "." + std::string(ext);
    writer->SetInputData(grid);
    writer->SetFileName(file_name.c_str());
    writer->SetDataModeToBinary();
    writer->Write();
    suanpan_debug("Plot is written to file \"{}\".\n", file_name);
}

void vtk_save_multiple(vtkNew<vtkMultiBlockDataSet>&& root, std::string file_name) {
    // NOLINT(performance-unnecessary-value-param)
    const vtkNew<vtkXMLMultiBlockDataWriter> writer;
    if(const auto ext = writer->GetDefaultFileExtension(); !file_name.ends_with(ext)) file_name += "." + std::string(ext);
    writer->SetInputData(root);
    writer->SetFileName(file_name.c_str());
    writer->SetDataModeToBinary();
    writer->Write();
    suanpan_debug("Plot is written to file \"{}\".\n", file_name);
}

int vtk_parser(const shared_ptr<DomainBase>& domain, std::istringstream& command) {
    auto plot_info = vtk_process(command);
    if(!plot_info.on_deformed) plot_info.scale = 0.;

    const auto L = to_token(to_category(plot_info.type));

    const auto func = OutputType::U == L || OutputType::V == L || OutputType::A == L || OutputType::RF == L || OutputType::DF == L || OutputType::IF == L ? vtk_plot_node_quantity : vtk_plot_element_quantity;

#ifdef SUANPAN_WIN
    domain->insert(std::async(std::launch::async, func, std::cref(domain), plot_info));
#else
    func(domain, plot_info);
#endif

    return SUANPAN_SUCCESS;
}

void vtk_plot_node_quantity(const shared_ptr<DomainBase>& domain, vtkInfo config) {
    auto& t_node_pool = domain->get_node_pool();
    auto& t_element_pool = domain->get_element_pool();

    config.title_name = "Plotting Nodal Quantity " + std::string(to_name(config.type));

    auto max_node = static_cast<unsigned>(t_node_pool.size());
    for(const auto& t_node : t_node_pool) max_node = std::max(max_node, t_node->get_tag());

    const vtkNew<vtkDoubleArray> data;
    const vtkNew<vtkPoints> node;

    data->SetNumberOfComponents(6);
    data->SetNumberOfTuples(++max_node);
    data->SetName(to_category(config.type).c_str());
    data->SetComponentName(0, "1");
    data->SetComponentName(1, "2");
    data->SetComponentName(2, "3");
    data->SetComponentName(3, "4");
    data->SetComponentName(4, "5");
    data->SetComponentName(5, "6");
    node->SetNumberOfPoints(max_node);

    for(auto I = 0u; I < max_node; ++I) data->SetTuple6(I, 0., 0., 0., 0., 0., 0.);

    vtkNew<vtkUnstructuredGrid> grid;
    grid->Allocate(static_cast<vtkIdType>(t_element_pool.size()));
    std::ranges::for_each(t_element_pool, [&](const shared_ptr<Element>& t_element) {
        const auto num_node = t_element->get_node_number();
        auto& encoding = t_element->get_node_encoding();

        if(const auto coor = t_element->GetDeformation(config.scale); coor.empty())
            for(auto I = 0u; I < num_node; ++I) node->SetPoint(encoding(I), 0., 0., 0.);
        else
            for(auto I = 0u; I < num_node; ++I) node->SetPoint(encoding(I), coor.colptr(I));

        t_element->GetData(data, to_token(to_category(config.type)));
        if(const auto t_cell = t_element->GetCell(); t_cell) grid->InsertNextCell(t_cell->GetCellType(), t_cell->GetPointIds());
    });

    grid->SetPoints(node);

    if(config.save_file) {
        grid->GetPointData()->SetScalars(data);
        domain->insert(std::async(std::launch::async, vtk_save_single, std::move(grid), config.file_name));
    }
    else {
        const vtkNew<vtkDoubleArray> sub_data;

        sub_data->SetNumberOfTuples(data->GetNumberOfTuples());
        sub_data->CopyComponent(0, data, to_index(config.type));

        grid->GetPointData()->SetScalars(sub_data);

        vtk_setup(grid, config);
    }
}

void vtk_plot_element_quantity_single(const shared_ptr<DomainBase>& domain, vtkInfo config) {
    auto& t_node_pool = domain->get_node_pool();
    auto& t_element_pool = domain->get_element_pool();

    config.title_name = "Plotting Element Quantity " + std::string(to_name(config.type));

    auto max_node = static_cast<unsigned>(t_node_pool.size());
    for(const auto& t_node : t_node_pool) max_node = std::max(max_node, t_node->get_tag());

    const vtkNew<vtkDoubleArray> data;
    const vtkNew<vtkPoints> node;

    data->SetNumberOfComponents(6);
    data->SetNumberOfTuples(++max_node);
    data->SetName(to_category(config.type).c_str());
    data->SetComponentName(0, "1");
    data->SetComponentName(1, "2");
    data->SetComponentName(2, "3");
    data->SetComponentName(3, "4");
    data->SetComponentName(4, "5");
    data->SetComponentName(5, "6");
    node->SetNumberOfPoints(max_node);

    for(unsigned I = 0; I < max_node; ++I) {
        node->SetPoint(I, 0., 0., 0.);
        data->SetTuple6(I, 0., 0., 0., 0., 0., 0.);
    }

    mat tensor(6, max_node, fill::zeros);
    vec counter(max_node, fill::zeros);

    vtkNew<vtkUnstructuredGrid> grid;
    grid->Allocate(static_cast<vtkIdType>(t_element_pool.size()));
    std::ranges::for_each(t_element_pool, [&](const shared_ptr<Element>& t_element) {
        if(-1 != config.material_type)
            if(const auto& t_tag = t_element->get_material_tag(); t_tag.empty() || static_cast<uword>(config.material_type) != t_tag(0)) return;

        const auto num_node = t_element->get_node_number();
        auto& encoding = t_element->get_node_encoding();

        if(const auto coor = t_element->GetDeformation(config.scale); coor.empty())
            for(auto I = 0u; I < num_node; ++I) node->SetPoint(encoding(I), 0., 0., 0.);
        else
            for(auto I = 0u; I < num_node; ++I) node->SetPoint(encoding(I), coor.colptr(I));

        counter(encoding) += 1.;
        if(const auto t_data = t_element->GetData(to_token(to_category(config.type))); !t_data.empty()) tensor.cols(encoding) += t_data;
        if(const auto t_cell = t_element->GetCell(); t_cell) grid->InsertNextCell(t_cell->GetCellType(), t_cell->GetPointIds());
    });

    for(auto I = 0u; I < max_node; ++I)
        if(0. != counter(I)) {
            tensor.col(I) /= counter(I);
            data->SetTuple(I, tensor.colptr(I));
        }

    grid->SetPoints(node);

    if(config.save_file) {
        grid->GetPointData()->SetScalars(data);
        domain->insert(std::async(std::launch::async, vtk_save_single, std::move(grid), config.file_name));
    }
    else {
        const vtkNew<vtkDoubleArray> sub_data;

        sub_data->SetNumberOfTuples(data->GetNumberOfTuples());
        sub_data->CopyComponent(0, data, to_index(config.type));

        grid->GetPointData()->SetScalars(sub_data);

        vtk_setup(grid, config);
    }
}

void vtk_plot_element_quantity_multiple(const shared_ptr<DomainBase>& domain, vtkInfo config) {
    auto& t_node_pool = domain->get_node_pool();
    auto max_node = static_cast<unsigned>(t_node_pool.size());
    for(const auto& t_node : t_node_pool) max_node = std::max(max_node, t_node->get_tag());

    const vtkNew<vtkPoints> global_node;
    global_node->SetNumberOfPoints(++max_node);

    config.title_name = "Plotting Element Quantity " + std::string(to_name(config.type));

    vtkNew<vtkMultiBlockDataSet> root;

    const auto category = to_category(config.type);
    const auto actual_type = to_token(category);

    auto counter{0u};

    for(const auto& t_element : domain->get_element_pool()) {
        if(-1 != config.material_type)
            if(const auto& t_tag = t_element->get_material_tag(); t_tag.empty() || static_cast<uword>(config.material_type) != t_tag(0)) continue;

        const auto num_node = t_element->get_node_number();
        auto encoding = t_element->get_node_encoding();

        if(const auto coor = t_element->GetDeformation(config.scale); coor.empty())
            for(auto I = 0u; I < num_node; ++I) global_node->SetPoint(encoding(I), 0., 0., 0.);
        else
            for(auto I = 0u; I < num_node; ++I) global_node->SetPoint(encoding(I), coor.colptr(I));

        const vtkNew<vtkPoints> node;
        node->SetNumberOfPoints(num_node);

        for(auto I = 0u; I < num_node; ++I) {
            node->SetPoint(I, global_node->GetPoint(static_cast<vtkIdType>(encoding(I))));
            encoding(I) = I;
        }

        const vtkNew<vtkDoubleArray> data;
        data->SetNumberOfComponents(6);
        data->SetNumberOfTuples(num_node);
        data->SetName(category.c_str());
        data->SetComponentName(0, "1");
        data->SetComponentName(1, "2");
        data->SetComponentName(2, "3");
        data->SetComponentName(3, "4");
        data->SetComponentName(4, "5");
        data->SetComponentName(5, "6");

        auto result = t_element->GetData(actual_type);
        if(result.empty()) result.resize(6, num_node);
        for(auto I = 0u; I < num_node; ++I) data->SetTuple(I, result.colptr(I));

        if(const auto cell = t_element->GetCell(encoding); cell) {
            const vtkNew<vtkUnstructuredGrid> grid;
            grid->Allocate(1);
            grid->InsertNextCell(cell->GetCellType(), cell->GetPointIds());
            grid->SetPoints(node);
            grid->GetPointData()->SetScalars(data);
            root->SetBlock(counter++, grid);
        }
    }

    if(config.save_file) domain->insert(std::async(std::launch::async, vtk_save_multiple, std::move(root), config.file_name));
}

class vtkBlock {
    const double scale;
    const OutputType field;
    const std::string name;

    const vtkNew<vtkPoints> grid_node;

public:
    const vtkNew<vtkUnstructuredGrid> grid;

    vtkBlock(const vtkIdType p, const double s, const OutputType f, const std::string_view n)
        : scale(s)
        , field(f)
        , name(n) { grid_node->SetNumberOfPoints(p); }

    auto process(const shared_ptr<Element>& element) const {
        const auto num_node = element->get_node_number();
        auto encoding = element->get_node_encoding();

        if(const auto coor = element->GetDeformation(scale); coor.empty())
            for(auto I = 0u; I < num_node; ++I) grid_node->SetPoint(encoding(I), 0., 0., 0.);
        else
            for(auto I = 0u; I < num_node; ++I) grid_node->SetPoint(encoding(I), coor.colptr(I));

        const vtkNew<vtkPoints> node;
        node->SetNumberOfPoints(num_node);

        for(auto I = 0u; I < num_node; ++I) {
            node->SetPoint(I, grid_node->GetPoint(static_cast<vtkIdType>(encoding(I))));
            encoding(I) = I;
        }

        const vtkNew<vtkDoubleArray> data;
        data->SetNumberOfComponents(6);
        data->SetNumberOfTuples(num_node);
        data->SetName(name.c_str());
        data->SetComponentName(0, "1");
        data->SetComponentName(1, "2");
        data->SetComponentName(2, "3");
        data->SetComponentName(3, "4");
        data->SetComponentName(4, "5");
        data->SetComponentName(5, "6");

        auto result = element->GetData(field);
        if(result.empty()) result.resize(6, num_node);
        for(auto I = 0u; I < num_node; ++I) data->SetTuple(I, result.colptr(I));

        if(const auto cell = element->GetCell(encoding); cell) {
            grid->InsertNextCell(cell->GetCellType(), cell->GetPointIds());
            grid->SetPoints(node);
            grid->GetPointData()->SetScalars(data);
        }
    }
};

void vtk_plot_element_quantity(const shared_ptr<DomainBase>& domain, vtkInfo config) {
    const auto func = config.average ? vtk_plot_element_quantity_single : vtk_plot_element_quantity_multiple;

    return func(domain, std::move(config));
}

#endif
